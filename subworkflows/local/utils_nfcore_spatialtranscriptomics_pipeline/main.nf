//
// Subworkflow with functionality specific to the nf-core/spatialtranscriptomics pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFVALIDATION_PLUGIN } from '../../nf-core/utils_nfvalidation_plugin'
include { paramsSummaryMap          } from 'plugin/nf-validation'
include { fromSamplesheet           } from 'plugin/nf-validation'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { dashedLine                } from '../../nf-core/utils_nfcore_pipeline'
include { nfCoreLogo                } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { workflowCitation          } from '../../nf-core/utils_nfcore_pipeline'

/*
========================================================================================
    SUBWORKFLOW TO INITIALISE PIPELINE
========================================================================================
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    help              // boolean: Display help text
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    pre_help_text = nfCoreLogo(monochrome_logs)
    post_help_text = '\n' + workflowCitation() + '\n' + dashedLine(monochrome_logs)
    def String workflow_command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"
    UTILS_NFVALIDATION_PLUGIN (
        help,
        workflow_command,
        pre_help_text,
        post_help_text,
        validate_params,
        "nextflow_schema.json"
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    emit:
    versions    = ch_versions
}

/*
========================================================================================
    SUBWORKFLOW FOR PIPELINE COMPLETION
========================================================================================
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:

    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(summary_params, email, email_on_fail, plaintext_email, outdir, monochrome_logs, multiqc_report.toList())
        }

        completionSummary(monochrome_logs)

        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

//
// Generate methods description for MultiQC
//
def toolCitationText() {

    def citation_text = [
            "Tools used in the workflow included:",
            "AnnData (Virshup et al. 2021),",
            "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016),",
            "Quarto (Allaire et al. 2022),",
            "Scanpy (Wolf et al. 2018),",
            "Space Ranger (10x Genomics)",
            "SpatialData (Marconato et al. 2023) and",
            "Squidpy (Palla et al. 2022)"
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {

    def reference_text = [
        '<li>Di Tommaso, P., Chatzou, M., Floden, E. W., Barja, P. P., Palumbo, E., & Notredame, C. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35(4), 316-319. doi: <a href="https://doi.org/10.1038/nbt.3820">10.1038/nbt.3820</a></li>',
        '<li>Ewels, P. A., Peltzer, A., Fillinger, S., Patel, H., Alneberg, J., Wilm, A., Garcia, M. U., Di Tommaso, P., & Nahnsen, S. (2020). The nf-core framework for community-curated bioinformatics pipelines. Nature Biotechnology, 38(3), 276-278. doi: <a href="https://doi.org/10.1038/s41587-020-0439-x">10.1038/s41587-020-0439-x</a></li>',
        '<li>Grüning, B., Dale, R., Sjödin, A., Chapman, B. A., Rowe, J., Tomkins-Tinch, C. H., Valieris, R., Köster, J., & Bioconda Team. (2018). Bioconda: sustainable and comprehensive software distribution for the life sciences. Nature Methods, 15(7), 475–476. doi: <a href="https://doi.org/10.1038/s41592-018-0046-7">10.1038/s41592-018-0046-7</a></li>',
        '<li>da Veiga Leprevost, F., Grüning, B. A., Alves Aflitos, S., Röst, H. L., Uszkoreit, J., Barsnes, H., Vaudel, M., Moreno, P., Gatto, L., Weber, J., Bai, M., Jimenez, R. C., Sachsenberg, T., Pfeuffer, J., Vera Alvarez, R., Griss, J., Nesvizhskii, A. I., & Perez-Riverol, Y. (2017). BioContainers: an open-source and community-driven framework for software standardization. Bioinformatics (Oxford, England), 33(16), 2580–2582. doi: <a href="https://doi.org/10.1093/bioinformatics/btx192">10.1093/bioinformatics/btx192</a></li>',
        '<li> Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA. bioRxiv 2021.12.16.473007. doi: <a href="https://doi.org/10.1101/2021.12.16.473007">10.1101/2021.12.16.473007</a></li>',
        '<li>Andrews S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data [Online]: <a href="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/">bioinformatics.babraham.ak.uk/project/fastqc</a></li>',
        '<li>Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016 Oct 1;32(19):3047-8. Epub 2016 Jun 16. PubMed PMID: 27312411; PubMed Central PMCID: PMC5039924. doi: <a href="https://doi.org/10.1093/bioinformatics/btw354">10.1093/bioinformatics/btw354</a></li>',
        '<li>Allaire J, Teague C, Scheidegger C, Xie Y, Dervieux C. Quarto (2022). doi: <a href="https://doi.org/10.5281/zenodo.5960048">10.5281/zenodo.5960048</a></li>',
        '<li>Wolf F, Angerer P, Theis F. SCANPY: large-scale single-cell gene expression data analysis. Genome Biol 19, 15 (2018). doi: <a href="https://doi.org/10.1186/s13059-017-1382-0">10.1186/s13059-017-1382-0</a></li>',
        '<li>10x Genomics Space Ranger 2.1.0 [Online]: <a href="https://www.10xgenomics.com/support/software/space-ranger">10xgenomics.com/support/software/space-ranger</a></li>',
        '<li>Marconato L, Palla G, Yamauchi K, Virshup I, Heidari E, Treis T, Toth M, Shrestha R, Vöhringer H, Huber W, Gerstung M, Moore J, Theis F, Stegle O. SpatialData: an open and universal data framework for spatial omics. bioRxiv 2023.05.05.539647; doi:<a href="https://doi.org/10.1101/2023.05.05.539647"> 10.1101/2023.05.05.539647</a></li>',
        '<li>Palla G, Spitzer H, Klein M et al. Squidpy: a scalable framework for spatial omics analysis. Nat Methods 19, 171–178 (2022). doi: <a href="https://doi.org/10.1038/s41592-021-01358-2">10.1038/s41592-021-01358-2</a></li>',
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
