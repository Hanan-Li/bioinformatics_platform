# This module takes in a workflow json file from the frontend and generates an R target config file that defines the pipeline.
import copy
import jinja2

# Global Dictionary.
NAME_TO_FUNCTION_DICT = {
    "read_seurat": "read_data",
    "get_genome_pattern": "append_genome_pattern",
    "v_plot": "violin_plot",
    "save_data": "persist_plot"
}
NAME_TO_FILE_STATUS_DICT = {
    "read_seurat": False,
    "save_data": True
}

# Gets the Relevant Jinja Template.
def getJinjaTemplate(template_file, template_directory):
    templateLoader = jinja2.FileSystemLoader(searchpath=template_directory)
    templateEnv = jinja2.Environment(loader=templateLoader)
    return templateEnv.get_template(template_file)

# Create parameters to pass into jinja to create DAG pipeline.
def createTargetFileParametersFromJson(dag):
    nodes = []
    for node in dag:
        processed_node = copy.deepcopy(node)
        processed_node["function"] = NAME_TO_FUNCTION_DICT[processed_node["name"]]
        processed_node["file"] = NAME_TO_FILE_STATUS_DICT.get(
            processed_node["name"], False)
        nodes.append(processed_node)
    return nodes

# Renders the R target config file.
def transformJsonAndCreateTargetsPipelineConfig(template_file, template_directory, params, outfile):
    template = getJinjaTemplate(template_file, template_directory)
    parameters = createTargetFileParametersFromJson(params)
    template.stream(nodes=parameters).dump(outfile)