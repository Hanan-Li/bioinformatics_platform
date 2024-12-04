import jwt

from datetime import datetime, timedelta, timezone
from functools import wraps
from flask import (
    request, jsonify, current_app
)

def token_required(f):
    @wraps(f)
    def decorated(*args, **kwargs):
        token = None
        # jwt is passed in the request header
        if 'x-access-token' in request.headers:
            token = request.headers['x-access-token']
        # return 401 if token is not passed
        if not token:
            return jsonify({'status': 'failure', 'error': 'Missing Token'}), 401

        try:
            # decoding the payload to fetch the stored details
            data = jwt.decode(token, current_app.config['SECRET_KEY'], algorithms="HS256")
            current_user = data['username']
            if datetime.now(timezone.utc).timestamp() >= data['exp']:
                return jsonify({'status': 'failure', 'error': 'Token Expired'}), 401
        except:
            return jsonify({'status': 'failure', 'error': 'Invalid Token'}), 401
        # returns the current logged in users context to the routes
        return f(current_user, *args, **kwargs)

    return decorated