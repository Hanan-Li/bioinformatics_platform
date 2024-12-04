import hashlib
import uuid
import jwt

from functools import wraps
from flask import (
    Blueprint, request, jsonify, Response, current_app
)
from datetime import datetime, timedelta, timezone
from api.db import get_db
from werkzeug.security import check_password_hash, generate_password_hash
bp = Blueprint('auth', __name__, url_prefix='/api/auth')



@bp.route('/register', methods=('POST',))
def register():
    if request.method == 'POST':
        data = request.get_json()
        email = data['email']
        username = data['username']
        password = data['password']
        db = get_db()
        error = None

        if not username:
            error = 'Username is required.'
        elif not password:
            error = 'Password is required.'
        elif not email:
            error = 'Email is required.'

        if error is None:
            try:
                db.execute(
                    "INSERT INTO user (email, username, password) VALUES (?, ?, ?)",
                    (email, username, sha512(password)),
                )
                db.commit()
            except db.IntegrityError:
                error = f"User {username} is already registered."
        if error is None:
            return jsonify({"status": "success", "username": username, "email": email, "error": ""}), 200, {'Content-Type': 'application/json'}
        else:
            return jsonify({"status": "failed", "username": "", "email": "", "error": error}), 500, {'Content-Type': 'application/json'}
    return jsonify({"status": "failed", "error": "method not allowed at endpoint"}), 405, {'Content-Type': 'application/json'}

# TODO: Return token instead of success status.


@bp.route('/login', methods=('POST',))
def login():
    data = request.get_json()
    username = data['username']
    print(data)
    input_password = data['password']
    db = get_db()
    db_row = db.execute(
        "SELECT * FROM user WHERE username = ?;", (username,)).fetchone()
    if db_row is not None:
        db_password = db_row['password']
        pass_list = db_password.split('$')
        hashed_input = pseudosha(pass_list[0], pass_list[1], input_password)
        if hashed_input == db_password:
            jwt_token = jwt.encode({'username': username, 'exp': datetime.now(
                tz=timezone.utc) + timedelta(minutes=30)}, current_app.config['SECRET_KEY'], algorithm="HS256")
            return jsonify({"status": "success", 'token': jwt_token}), 200, {'Content-Type': 'application/json'}
        else:
            return jsonify({"status": "failure", "error": "incorrect password"}), 401, {'Content-Type': 'application/json'}
    return jsonify({"status": "failure", "error": "user does not exist"}), 401, {'Content-Type': 'application/json'}


def pseudosha(algorithm, salt, password):
    """pseudosha."""
    hash_obj = hashlib.new(algorithm)
    password_salted = salt + password
    hash_obj.update(password_salted.encode('utf-8'))
    password_hash = hash_obj.hexdigest()
    password_db_string = "$".join([algorithm, salt, password_hash])
    return password_db_string


def sha512(password):
    """sha512."""
    algorithm = 'sha512'
    salt = uuid.uuid4().hex
    hash_obj = hashlib.new(algorithm)
    password_salted = salt + password
    hash_obj.update(password_salted.encode('utf-8'))
    password_hash = hash_obj.hexdigest()
    password_db_string = "$".join([algorithm, salt, password_hash])
    return password_db_string
