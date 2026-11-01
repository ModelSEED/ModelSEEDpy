"""Unit tests for the JSON-RPC client's error handling.

Regression coverage for the 500-response branch. Different upstream
services return different error-body shapes (JSON-RPC servers should
return a dict, but some backends occasionally return a string or a
list under load). The client must translate all of them into a
`ServerError` rather than crashing with a `TypeError` from `**` on
a non-mapping.
"""

from __future__ import annotations

import pytest

from modelseedpy.core.rpcclient import RPCClient, ServerError


class _FakeResp:
    def __init__(self, status_code: int, headers: dict, body):
        self.status_code = status_code
        self.headers = headers
        self.text = body if isinstance(body, str) else ""
        self._body = body

    @property
    def encoding(self):
        return "utf-8"

    @encoding.setter
    def encoding(self, _value):
        pass

    @property
    def ok(self):
        return self.status_code < 400

    def json(self):
        if isinstance(self._body, (dict, list)):
            return self._body
        import json as _json
        return _json.loads(self._body)

    def raise_for_status(self):
        if not self.ok:
            raise RuntimeError(f"HTTP {self.status_code}")


def _client() -> RPCClient:
    return RPCClient("http://example.invalid/rpc")


def _patched(body_shape) -> _FakeResp:
    """Build the FakeResp for a 500 JSON response with the given error shape."""
    return _FakeResp(
        status_code=500,
        headers={"content-type": "application/json"},
        body={"error": body_shape},
    )


def test_500_with_dict_error_raises_serverererror_populated(monkeypatch):
    body = {"name": "JSONRPCError", "code": -32603, "message": "boom", "error": "trace-here"}
    monkeypatch.setattr(
        "modelseedpy.core.rpcclient._requests.post",
        lambda *a, **kw: _patched(body),
    )
    with pytest.raises(ServerError) as exc_info:
        _client().call("some.method", [])
    assert exc_info.value.name == "JSONRPCError"
    assert exc_info.value.code == -32603
    assert exc_info.value.message == "boom"


def test_500_with_string_error_does_not_typeerror(monkeypatch):
    """Regression: tutorial.theseed.org and some other backends return
    a plain string here under load. Pre-fix, this raised
    `TypeError: argument after ** must be a mapping, not str` and
    masked the real upstream failure. Post-fix, the string is
    surfaced as a ServerError message."""
    monkeypatch.setattr(
        "modelseedpy.core.rpcclient._requests.post",
        lambda *a, **kw: _patched("service temporarily unavailable"),
    )
    with pytest.raises(ServerError) as exc_info:
        _client().call("some.method", [])
    assert exc_info.value.name == "Unknown"
    assert exc_info.value.code == 0
    assert "service temporarily unavailable" in exc_info.value.message


def test_500_with_list_error_does_not_typeerror(monkeypatch):
    """Same as the string case but for a list body (also non-mapping)."""
    monkeypatch.setattr(
        "modelseedpy.core.rpcclient._requests.post",
        lambda *a, **kw: _patched(["a", "b"]),
    )
    with pytest.raises(ServerError) as exc_info:
        _client().call("some.method", [])
    assert exc_info.value.name == "Unknown"
    assert "['a', 'b']" in exc_info.value.message


def test_500_with_no_error_key_raises_unknown(monkeypatch):
    """Body is JSON but has no 'error' field. Pre-existing behavior:
    raise ServerError('Unknown', 0, ret.text). Locking it in as a
    guard against regression."""
    monkeypatch.setattr(
        "modelseedpy.core.rpcclient._requests.post",
        lambda *a, **kw: _FakeResp(500, {"content-type": "application/json"}, {"result": None}),
    )
    with pytest.raises(ServerError) as exc_info:
        _client().call("some.method", [])
    assert exc_info.value.name == "Unknown"


def test_500_with_non_json_body_raises_unknown(monkeypatch):
    """Content-Type isn't JSON. Pre-existing branch. Regression-locked."""
    monkeypatch.setattr(
        "modelseedpy.core.rpcclient._requests.post",
        lambda *a, **kw: _FakeResp(500, {"content-type": "text/html"}, "<h1>Bad Gateway</h1>"),
    )
    with pytest.raises(ServerError) as exc_info:
        _client().call("some.method", [])
    assert exc_info.value.name == "Unknown"
    assert "Bad Gateway" in exc_info.value.message
