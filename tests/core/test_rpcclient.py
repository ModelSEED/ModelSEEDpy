# -*- coding: utf-8 -*-
import json
import pytest
import modelseedpy.core.rpcclient
from modelseedpy.core.rpcclient import RPCClient, ServerError


class MockResponse:
    def __init__(self, body):
        self.status_code = 500
        self.ok = False
        self.headers = {"content-type": "application/json"}
        self.text = json.dumps(body)
        self._body = body

    def json(self):
        return self._body


def mock_post(body):
    return lambda *args, **kwargs: MockResponse(body)


def test_call_server_error(monkeypatch):
    error = {"name": "JSONRPCError", "code": -32603, "message": "bad input"}
    monkeypatch.setattr(
        modelseedpy.core.rpcclient._requests, "post", mock_post({"error": error})
    )
    with pytest.raises(ServerError) as e:
        RPCClient("http://localhost").call("Service.method", [])
    assert e.value.name == "JSONRPCError"
    assert e.value.code == -32603
    assert e.value.message == "bad input"


@pytest.mark.parametrize(
    "body, name, code, message",
    [
        # JSON-RPC 2.0 error object, which has no "name"
        (
            {"error": {"code": -32601, "message": "Method not found"}},
            "Unknown",
            -32601,
            "Method not found",
        ),
        # error object with a key ServerError does not take
        (
            {"error": {"name": "Err", "code": 1, "message": "m", "stack": "..."}},
            "Err",
            1,
            "m",
        ),
        # error that is not an object
        ({"error": "service unavailable"}, "Unknown", 0, "service unavailable"),
        # body that is not an object
        ("internal error", "Unknown", 0, "internal error"),
    ],
)
def test_call_server_error_other_shapes(monkeypatch, body, name, code, message):
    """
    A 500 response must raise ServerError whatever the shape of its error body
    """
    monkeypatch.setattr(modelseedpy.core.rpcclient._requests, "post", mock_post(body))
    with pytest.raises(ServerError) as e:
        RPCClient("http://localhost").call("Service.method", [])
    assert e.value.name == name
    assert e.value.code == code
    assert message in e.value.message
