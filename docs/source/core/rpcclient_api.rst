rpcclient
---------------------------

+++++++++++++++++++++
ServerError()
+++++++++++++++++++++

This class constructs a descriptive string of a Server Error:

.. code-block:: python

 from modelseedpy.core import ServerError
 servErr = ServerError(name, code, message, data=None, error=None)

- *name* & *message* ``str``: The name and message of the error.
- *code* ``int``: The numerical code number of the error.
- *data* & *error* ``str``: The data and error of the server error, respectively, where only one of these will be included when the string is constructed.

----------------
__str__()
----------------

**returns** *error_string* ``str``: A descriptive string of the Server Error that consists of the content that was provided in the initiation of the class.

-------------------------------------
Accessible content
-------------------------------------

The *name*, *message*, *code*, and *data* or *error* parameters that comprise the *error_string* from the ``__str__ function are accessible from the class.

+++++++++++++++++++++
RPCClient()
+++++++++++++++++++++

This class offers a suite of static method functions that assist users in editing and expanding COBRA models:

.. code-block:: python

 from modelseedpy.core import RPCClient
 rpcCli = RPCClient(url,token=None,version="1.0",timeout=30 * 60,trust_all_ssl_certificates=False)

- *url* & *token* ``str``: The URL that will be parsed via requests, and the token that authorizes access to the URL site.
- *version* ``str``: The version of the URL that will be parsed.
- *timeout* ``float``: The limit of seconds at which requests post will terminate.
- *trust_all_ssl_certificates* ``bool``: specifies whether the requests posts will be verified.

-------------
call()
-------------

**returns** *error_string* ``str``: A descriptive string of the Server Error that consists of the content that was provided in the initiation of the class.

.. code-block:: python

 resp_result = rpcCli.call(method, params, token=None)

- *method* ``str`` & *params* ``dict``: Components that are passed as data to the URL request.
- *token* ``str``: The token that authorizes access to the URL site.

**returns** *resp_result* ``str``: The requests result where it exists, otherwise ``None``.

-------------------------------------
Accessible content
-------------------------------------

The *url*, *version*, *token*, *timeout*, and *trust_all_ssl_certificates* parameters are accessible from the class.
