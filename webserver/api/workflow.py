import functools
import os
import json
import random
import string

from flask import (
    Blueprint, flash, g, redirect, render_template, request, session, url_for, jsonify, current_app
)
from werkzeug.security import check_password_hash, generate_password_hash
from api.common import token_required
from api.db import get_db

import requests

bp = Blueprint('workflow', __name__, url_prefix='/api/workflow')

@bp.route('/create', methods=('GET', 'POST'))
@token_required
def create_workflow(username):
    if request.method == 'GET':
        print(username)
        return jsonify({"test": "test"})
    return ''


@bp.route('/workflow_status', methods=('GET',))
def workflow_status():
    return ''

@bp.route('/get', methods=('GET',))
@token_required
def get_workflows(username):
    print(username)
    dir_path = os.path.dirname(os.path.realpath(__file__))
    pipeline_defs = os.path.join(
        dir_path, current_app.config['PIPELINE_DEFINITIONS_DIR'])
    pipeline_filenames = [os.path.join(pipeline_defs, f) for f in os.listdir(
        pipeline_defs) if os.path.isfile(os.path.join(pipeline_defs, f)) and f != 'README.md']
    print(pipeline_filenames)
    pipelines = []
    for pipeline_json in pipeline_filenames:
        with open(pipeline_json) as f:
            content = f.read()
            pipelines.append(json.loads(content))
    print(pipelines)
    return jsonify(pipelines), 200


@bp.route('/instantiate', methods=('POST',))
@token_required
def instantiate_workflows(username):
    data = request.get_json()
    # Hardcode static value for now
    runName = ''.join(random.choices(
        string.ascii_letters, k=12))
    payload = f'{{"launch":{{"id":"LU58zQrxnDvVjwMV11D7B","computeEnvId":"6ZVpe9ADcH2Kvrx6qlrdFF","pipeline":"https://github.com/nextflow-io/rnaseq-nf","pipelineId":74554027315001,"revision":null,"workDir":"gs://seqeratest","configProfiles":["standard"],"configText":null,"towerConfig":null,"paramsText":"{{\\"outdir\\": \\"results\\",  \\"reads\\": \\"${{projectDir}}/data/ggal/ggal_gut_{{1,2}}.fq\\",  \\"multiqc\\": \\"${{projectDir}}/multiqc\\",  \\"transcriptome\\": \\"${{projectDir}}/data/ggal/ggal_1_48850000_49020000.Ggal71.500bpflank.fa\\"}}","preRunScript":null,"postRunScript":null,"pullLatest":false,"stubRun":false,"mainScript":null,"entryName":null,"schemaName":null,"workspaceSecrets":[],"userSecrets":[],"runName":"{runName}","labelIds":[],"optimizationId":null}}}}'
    response = callSeqera("workflow/launch/", payload, "POST").json()
    print(data)
    print(response)
    db = get_db()
    db_row = db.execute(
        "SELECT * FROM user WHERE username = ?;", (username,)).fetchone()
    if db_row is not None:
        user_id = db_row['id']
        db.execute(
            "INSERT INTO instantiation (user_id, workflow_id, project_id, seqera_workflow_id) VALUES (?, ?, ?, ?)",
            (user_id, data["workflow_id"], 0, response["workflowId"]),
        )
        db.commit()
        return jsonify({"seqera_workflow_id": response["workflowId"], "status": "success"}), 200
    return jsonify({"test": "test"}), 200


def callSeqera(api_path, payload, method):
    full_api_path = "https://api.cloud.seqera.io/" + api_path
    print(full_api_path)
    headers = {"Authorization": "Bearer " + current_app.config['SEQERA_TOKEN']}
    if method == "GET":
        return requests.get(full_api_path, headers=headers)
    else:
        return requests.post(full_api_path, json=json.loads(payload), headers=headers)

