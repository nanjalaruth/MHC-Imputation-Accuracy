#!/usr/bin/env nextflow

/*
========================================================================================
=                                 h3achipimputation                                    =
========================================================================================
 h3achipimputation imputation pipeline.
----------------------------------------------------------------------------------------
 @Homepage / @Documentation
 
https://github.com/nanjalaruth/MHC-Imputation-Accuracy
----------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------
================================================================================
=                           C O N F I G U R A T I O N                          =
================================================================================
*/

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

// Configurable variables
params.name = false
params.email = false
params.plaintext_email = false

output_docs = file("$baseDir/docs/output.md")

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
    custom_runName = workflow.runName
}

// Header log info
log.info """
=======================================================
HLAimputation v${params.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']    = 'HLAimputation'
summary['Pipeline version'] = params.version
summary['Run Name']         = custom_runName ?: workflow.runName
summary['Target datasets']  = params.target_datasets.values().join(', ')
summary['Reference panels']  = params.ref_panels.keySet().join(', ')
summary['Max Memory']       = params.max_memory
summary['Max CPUs']         = params.max_cpus
summary['Max Time']         = params.max_time
summary['Output dir']       = params.outDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['Current path']     = "$PWD"
summary['Git info']         = "${workflow.repository} - ${workflow.revision} [${workflow.commitId}]"
summary['Command line']     = workflow.commandLine
if(workflow.containerEngine) {
    summary['Container Engine'] = workflow.containerEngine
    summary['Container'] = workflow.container
    summary['Current home'] = "$HOME"
    summary['Current user'] = "$USER"
    summary['Current path'] = "$PWD"
    summary['Working dir'] = workflow.workDir
    summary['Output dir'] = params.outDir
    summary['Script dir'] = workflow.projectDir
    summary['Config Profile'] = workflow.profile
}

if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


def create_workflow_summary(summary) {

    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'h3achipimputation-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'h3achipimputation Workflow Summary'
    section_href: 'https://github.com/h3abionet/chipimputation'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

    return yaml_file
}

