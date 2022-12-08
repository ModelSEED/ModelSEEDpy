# -*- coding: utf-8 -*-

from __future__ import absolute_import

import json as _json
import requests as _requests
import random as _random


class _JSONObjectEncoder(_json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, set):
            return list(obj)
        if isinstance(obj, frozenset):
            return list(obj)
        return _json.JSONEncoder.default(self, obj)


class ServerError(Exception):
    def __init__(self, name, code, message=None, data=None, error=None):
        super(Exception, self).__init__(message)
        self.name = name
        self.code = code
        self.message = message or ""
        self.data = data or error or ""
        # data = JSON RPC 2.0, error = 1.1

    def __str__(self):
        return (
            self.name + ": " + str(self.code) + ". " + self.message + "\n" + self.data
        )


class RPCClient:
    def __init__(
        self,
        url,
        token=None,
        version="1.0",
        timeout=30 * 60,
        trust_all_ssl_certificates=False,
    ):
        self.url = url
        self.token = token
        self.version = version
        self.timeout = timeout
        self.trust_all_ssl_certificates = trust_all_ssl_certificates

    def call(self, method, params, token=None):
        headers = {}
        if token:
            headers["AUTHORIZATION"] = token
        elif self.token:
            headers["AUTHORIZATION"] = self.token
        arg_hash = {
            "method": method,
            "params": params,
            "version": self.version,
            "id": str(_random.random())[2:],
            "context": {},
        }
        body = _json.dumps(arg_hash, cls=_JSONObjectEncoder)
        ret = _requests.post(
            self.url,
            data=body,
            headers=headers,
            timeout=self.timeout,
            verify=not self.trust_all_ssl_certificates,
        )
        ret.encoding = "utf-8"
        if ret.status_code == 500:
            if ret.headers.get("content-type") == "application/json":
                err = ret.json()
                if "error" in err:
                    raise ServerError(**err["error"])
                else:
                    raise ServerError("Unknown", 0, ret.text)
            else:
                raise ServerError("Unknown", 0, ret.text)
        if not ret.ok:
            ret.raise_for_status()
        resp = ret.json()
        if "result" not in resp:
            raise ServerError("Unknown", 0, "An unknown server error occurred")
        if not resp["result"]:
            return None
        return resp["result"]
