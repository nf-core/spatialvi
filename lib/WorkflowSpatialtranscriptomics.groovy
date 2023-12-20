//
// This file holds several functions specific to the workflow/spatialtranscriptomics.nf in the nf-core/spatialtranscriptomics pipeline
//

import nextflow.Nextflow
import groovy.text.SimpleTemplateEngine

class WorkflowSpatialtranscriptomics {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {

        // if (!params.fasta) {
        //     Nextflow.error "Genome fasta file not specified with e.g. '--fasta genome.fa' or via a detectable config file."
        // }
    }

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }

    //
    // Generate methods description for MultiQC
    //

    public static String toolCitationText(params) {

        def citation_text = [
                "Tools used in the workflow included:",
                "AnnData (Virshup et al. 2021),",
                "FastQC (Andrews 2010),",
                "MultiQC (Ewels et al. 2016),",
                "Quarto (Allaire et al. 2022),",
                "Scanpy (Wolf et al. 2018),",
                "Space Ranger (10x Genomics) and",
                "SpatialDE (Svensson et al. 2018)."
            ].join(' ').trim()

        return citation_text
    }

    public static String toolBibliographyText(params) {

        def reference_text = [
                "<li>Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA. bioRxiv 2021.12.16.473007; doi: <a href=https://doi.org/10.1101/2021.12.16.473007>10.1101/2021.12.16.473007</a></li>",
                "<li>Andrews S, (2010) FastQC, URL: <a href=https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>bioinformatics.babraham.ac.uk</a>.</li>",
                "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: <a href=https://doi.org/10.1093/bioinformatics/btw354>10.1093/bioinformatics/btw354</a></li>",
                "<li>Allaire J, Teague C, Scheidegger C, Xie Y, Dervieux C. Quarto (2022). doi: <a href=https://doi.org/10.5281/zenodo.5960048>10.5281/zenodo.5960048</a></li>",
                "<li>Wolf F, Angerer P, Theis F. SCANPY: large-scale single-cell gene expression data analysis. Genome Biol 19, 15 (2018). doi: <a href=https://doi.org/10.1186/s13059-017-1382-0>10.1186/s13059-017-1382-0</a></li>",
                "<li>10x Genomics Space Ranger 2.1.0, URL: <a href=https://www.10xgenomics.com/support/software/space-ranger>10xgenomics.com/support/software/space-ranger</a></li>",
                "<li>Svensson V, Teichmann S, Stegle O. SpatialDE: identification of spatially variable genes. Nat Methods 15, 343–346 (2018). doi: <a href=https://doi.org/10.1038/nmeth.4636>10.1038/nmeth.4636</a></li>",
            ].join(' ').trim()

        return reference_text
    }

    public static String methodsDescriptionText(run_workflow, mqc_methods_yaml, params) {
        // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
        def meta = [:]
        meta.workflow = run_workflow.toMap()
        meta["manifest_map"] = run_workflow.manifest.toMap()

        // Pipeline DOI
        meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
        meta["nodoi_text"] = meta.manifest_map.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

        // Tool references
        meta["tool_citations"] = ""
        meta["tool_bibliography"] = ""
        meta["tool_citations"] = toolCitationText(params).replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
        meta["tool_bibliography"] = toolBibliographyText(params)


        def methods_text = mqc_methods_yaml.text

        def engine =  new SimpleTemplateEngine()
        def description_html = engine.createTemplate(methods_text).make(meta)

        return description_html
    }

}
