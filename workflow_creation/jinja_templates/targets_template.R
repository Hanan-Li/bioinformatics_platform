library(targets)
tar_source("functions/seurat_functions.R")
tar_option_set(packages = c("Seurat", "Matrix", "irlba", "RSpectra", "ggplot2", "sctransform", "irlba"), error="null")

list(
    {% for node in nodes %}
    tar_target({{ node.name }}, {{ node.function}}({% for param in node.params %}{{param}}{% if not loop.last %}, {% endif %}{% endfor %}){% if node.file %}, format = 'file'{% endif %}){% if not loop.last %}, {% endif %}
    {% endfor %}
)