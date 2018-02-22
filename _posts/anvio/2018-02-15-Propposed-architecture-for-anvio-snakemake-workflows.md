---
layout: post
title: Proposed architecture for anvio snakemake workflows
modified: 2018-02-15
excerpt: "Bringing the magic of anvio together with the wonder of snakemake."
comments: true
authors: [alon]
categories: [anvio]
---

{% include _toc.html %}

# Proposed architecture for anvi'o snakemake workflows

snakemake is a very useful and robust language for creating computational pipelines. Recently we have been using it extensivley and it has been wonderful.

In order to let others enjoy the wonder of anvi'o together with the wonder of snakemake, we embarked on an effort to make "easy-to-use" (well, not TOO easy, after all this is science) snakefiles for some of the commonly used workflows involving anvi'o.

In this post I wish to propose an architecture for how these types of workflows should be written. I have been developing workflows according to this architecture in the branch `workflows` of anvi'o. A standard architecture would allow us to automate the writing of new workflows. It would also make it easier to understand and change existing ones. Finally, it would make it easier to develope user friendly warning and error messages so that the user could enjoy all the flexibility of the workflows.

When developing these workflows I made an effort to allow the user to easily configure any configurable parameter that is accessible by the software that are utilized in the workflow. This was not always possible to do (under reasonable effort), as some constraints have to be made for the workflow to, well, work and flow.

## Keeping the workflows up-to-date with anvi'o versions

By placing the workflows in the anvi'o repository, and updating them whenever a change is required due to a change in the anvi'o programs, we can guarantee that the user would always have workflows that are compatible with the version of anvi'o.

## The `anvio.workflows` package

The `anvio.workflows` package contains modules to help implement the workflows. The main purpose of these modules is to deal with configurable parameters of the rules in the workflow. By dealing, I mean enabling easy access, and sanity checks. For convinience I always include `import anvio.workflows as w` at the head of each snakefile.

## The `anvio.workflows.workflowsops` module

Each workflow corresponds to a class in `anvio.workflows.workflowsops`. Each of these classes inherits from `WorkflowSuperClass`, and in addition it inherits from
any class corresponding to a workflow that is included in the workflow (i.e. if
there is an `include` statement in the workflow). For example the `pangenome.snake` has the following line:

```
include: w.get_path_to_workflows_dir() + "/generate_and_annotate_contigs_db.snake"`
```

Which means it includes the all the rules from  `generate_and_annotate_contigs_db.snake`, and accordingly the class `PangenomicsWorkflow` inherits from `ContigsDBWorkflow`.

An object of `WorkflowSuperClass` has the following attributes:

	rules - a list of all the rules in the workflow
	config - the config file dictionary
	acceptable_params_dict - a dictionary with the rule names as keys. Each dictionaty entry is a list of the acceptable patameters for a rule, i.e. the parameters that could be changed via the config file.
	dirs_dict - see below in "Definitions of directories in the workflow"
	general_params - other workflow related parameters that could be found in the config file. These are all the parameters that are not related to just one specific rule (hence, "general"), and as such, they will not be nested under a rule name in the config file.

	
This allows us to easily check if the config file contains parameters that are not used in the workflow (and hence raise an error and notify the user). In addition, each instance of `WorkflowSuperClass`, can have an empty config file genearted, which contains all the configurable parameters for the workflow (by running `self.save_empty_config_in_json_format()`. We can later rap this nicely with something like `anvi-gen-empty-config --workflow NAMEOFWORKFLOW --config-file-name CONFIGFILENAME`, and then the user would be able to simply delete the parameters they don't need, and define the parameters they want. If we are truly up to it, then we would make another module `get_default_config`, but for now it was just too much work.

For example:

The usage of these classes is as follows:

```python
from anvio.workflows.workflowsops import PangenomicsWorkflow
p = PangenomicsWorkflow(config)
p.init()
p.save_empty_config_in_json_format("empty-config-for-pangenomics.json")
```

### Definitions of directories in the workflow

In order to allow the user to define the names of directories through the config file, any directory that is expected to be created by the workflow must be defined in the `dirs_dict` in `anvio.workflows`, with a default name. For now this is what we have:

```python
dirs_dict = {"LOGS_DIR"     : "00_LOGS"         ,\
             "QC_DIR"       : "01_QC"           ,\
             "ASSEMBLY_DIR" : "02_ASSEMBLY"     ,\
             "CONTIGS_DIR"  : "03_CONTIGS"      ,\
             "MAPPING_DIR"  : "04_MAPPING"      ,\
             "PROFILE_DIR"  : "05_ANVIO_PROFILE",\
             "MERGE_DIR"    : "06_MERGED"       ,\
             "PAN_DIR"      : "07_PAN"          ,\
             "FASTA_DIR"    : "01_FASTA"        ,\
             "LOCI_DIR"     : "04_LOCI_FASTAS"   \
}
```

The user-defined names from the config file, are set when the `WorkflowSuperClass.init()` is ran (see `p.init()` in the example above).

### The number of threads for a rule

In order to use the functionality of `threads:` in snakemake, each rule has the following lines included in the definition:

```
    threads: w.T(config, 'RULENAME', 1)
    resources: nodes = w.T(config, 'RULENAME', 1)
```

## Alon's rules for Defining rules

First of all, the default target rule (i.e. the first rule in the snakefile) is NOT called "all". Simply because if it were then when we use `includes:` (see it [here](http://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#includes) in the snakemake docummentation), then we would risk having multiple rules with the same name (which is not allowed).

This is the template I use for rules:

```
rule NEWRULE:
    version: 1.0
    log: dirs_dict["LOGS_DIR"] + "/{sample}-NEWRULE.log"
    input:
    output:
    params:
    threads: w.T(config, 'NEWRULE', 1)
    resources: nodes = w.T(config, 'NEWRULE', 1)
    shell: " >> {log} 2>&1"
```

I will use the rule for `anvi-pan-genome` from the pangenomic workflow as an example of how I write rules, and utilize the `anvio.workflows` package.

If a rule is basically wrapping a single program, then I set the name of the rule to be the same as the name of the program (except that any `-` has to be replaced by `_`). Below you can see the rule `anvi_pan_genome`:

```python
rule anvi_pan_genome:
    version: anvio.__pan__version__
    log: dirs_dict["LOGS_DIR"] + "/anvi_pan_genome.log"
    threads: w.T(config, "anvi_pan_genome", 20)
    resources: nodes = w.T(config, "anvi_pan_genome", 20)
    input: dirs_dict["PAN_DIR"] + "/" + project_name + "-GENOMES.db"
    params:
        output_dir = dirs_dict["PAN_DIR"],
        genome_names = w.B(config, "anvi_pan_genome", "--genome-names"),
        project_name = pan_project_name,
        skip_alignments = w.B(config, "anvi_pan_genome", "--skip-alignments"),
        align_with = w.B(config, "anvi_pan_genome", "--align-with"),
        exclude_partial_gene_calls = w.B(config, "anvi_pan_genome", "--exclude-partial-gene-calls"),
        use_ncbi_blast = w.B(config, "anvi_pan_genome", "--use-ncbi-blast"),
        minbit = w.B(config, "anvi_pan_genome", "--minbit"),
        mcl_inflation = w.B(config, "anvi_pan_genome", "--mcl-inflation"),
        min_occurrence = w.B(config, "anvi_pan_genome", "--min-occurrence"),
        min_percent_identity = w.B(config, "anvi_pan_genome", "--min-percent-identity"),
        sensitive = w.B(config, "anvi_pan_genome", "--sensitive"),
        description = w.B(config, "anvi_pan_genome", "--description"),
        overwrite_output_destinations = w.B(config, "anvi_pan_genome", "--overwrite-output-destinations"),
        skip_hierarchical_clustering = w.B(config, "anvi_pan_genome", "--skip-hierarchical-clustering"),
        enforce_hierarchical_clustering = w.B(config, "anvi_pan_genome", "--enforce-hierarchical-clustering"),
        distance = w.B(config, "anvi_pan_genome", "--distance"),
        linkage = w.B(config, "anvi_pan_genome", "--linkage")
    output: dirs_dict["PAN_DIR"] + "/" + pan_project_name + "-PAN.db"
    shell:
        """
            anvi-pan-genome -g {input} --num-threads {threads} -o {params.output_dir} {params.genome_names}\
            {params.skip_alignments} {params.align_with} {params.exclude_partial_gene_calls}\
            {params.use_ncbi_blast} {params.minbit} {params.mcl_inflation}\
            {params.min_occurrence} {params.min_percent_identity} {params.sensitive}\
            {params.project_name} {params.description} {params.overwrite_output_destinations}\
            {params.skip_hierarchical_clustering} {params.enforce_hierarchical_clustering}\
            {params.distance} {params.linkage}
        """
```

Ok so it looks like a lot, so let's dive in... In an effort to make things flexible I allow the user to change all parameters that are acceptable to `anvi-pan-genome`. In order to make things easy, I just named the parameters the same as the "double-dash" version of the arguments. This works for anvi'o most of the time as we make an effort to provide a comprehensive name to each argument.

To explain the module `w.B`, works slightly differently for arguments that are flags (such as `--use-ncbi-blast`) and arguments that accept an input (such as `--min-occurrence`). In both cases, if the argument is not defined in the config file, then it would return either an empty string or the default that was provided to `w.B`. For flags, it would simply return the argument (for example, `"--use-ncbi-blast"`), but for arguments such as `--min-occurrence` it would return `--min-occurrence VALUE` where `VALUE` is the value provided in the config file. That's why in the shell command it is enough to write `{params.min_occurrence}`, and there is no reason to write `--min-occurrence {params.min_occurrence}`.
