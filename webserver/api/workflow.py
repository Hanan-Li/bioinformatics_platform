import functools

from flask import (
    Blueprint, flash, g, redirect, render_template, request, session, url_for, jsonify
)
from werkzeug.security import check_password_hash, generate_password_hash
from api.common import token_required
from . import targets_pipeline_parser
bp = Blueprint('workflow', __name__, url_prefix='/api/workflow')

@bp.route('/create', methods=('GET', 'POST'))
@token_required
def create_workflow(username):
    if request.method == 'GET':
        print(username)
        return jsonify({"test": "test"})
    
    data = request.json
    targets_pipeline_parser.transformJsonAndCreateTargetsPipelineConfig("targets_template.R", "/home/hananli/Documents/bioinformatics_platform/workflow_creation/jinja_templates", data)
    return ''

# @bp.route('/workflow_status', methods= ('GET'))
# def workflow_status():
#     return ''