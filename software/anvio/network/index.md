---
title: The anvi'o universe
excerpt: "Inputs, outputs, programs, and concepts. All linked. Just like the way you like them."
layout: project-anvio-network
modified: 2018-06-23
categories: [anvio]
redirect_from: /nt/
image:
  feature: anvio-network.png
  display: false
---

<div id="infopanel">

<div id="logo">
    <a href="https://github.com/merenlab/anvio/releases" target="_blank">
        <img src="https://github.com/merenlab/anvio/raw/master/anvio/data/interactive/images/logo.png" style="width:222px; opacity: 0.8;" />
    </a>
</div>
<p id="description">The purpose of this page is to visualize <b>the interconnected nature of anvi'o concepts and programs</b>. Visit our <b><a href="/software/anvio/help/" target="_blank">help pages</a></b> for more information on anvi'o programs and concepts.
<br /><br />
Through its modular design principles, <b>anvi'o resembles LEGO, where the same simple blocks can build non-identical complex structures</b>. Interfaces between anvi'o tools are established by self-contained databases that are generated, modified, queried, visualized, and merged with other databases through atomic anvi'o programs. Through this modularity, anvi'o aims to empower microbiologists without imposing rigid workflows or adding substantial limitations on creativity and to help them explore large datasets through their unique and diverse approaches.
<br /><br />
You can zoom-in, zoom-out and move the network around. If you click on any item on this page, you will see how to acquire that item and what you can do with it: red lines show where is it coming from, and green lines show where can you go with it. The descriptions of each icon is below.
</p>

<div class="info">
<img src="/images/icons/PROGRAM.png" width="30" />
<p>An anvi'o program. Click info button next to program name to access to the help menu. A full list of anvi'o programs can be found <a href="/software/anvio/vignette/" target="_blank">here</a>.</p>
</div>

<div class="info">
<img src="/images/icons/DB.png" width="30" />
<p>An anvi'o database (special data storage files generated and accessed by anvi'o programs).</p>
</div>

<div class="info">
<img src="/images/icons/CONCEPT.png" width="30" />
<p>An anvi'o concept. Things that are meaningful within anvi'o, generated or used by anvi'o programs.</p>
</div>

<div class="info">
<img src="/images/icons/DISPLAY.png" width="30" />
<p>An anvi'o display for interactive data analysis and visualization.</p>
</div>

<div class="info">
<img src="/images/icons/TXT.png" width="30" />
<p>A TAB-delimited text file.</p>
</div>

<div class="info">
<img src="/images/icons/FASTA.png" width="30" />
<p>FASTA-formatted file for sequences in DNA or amino acid alphabet.</p>
</div>

<div class="info">
<img src="/images/icons/BAM.png" width="30" />
<p>Binary SAM, a common file format for widely used mapping software to report short read alignments.</p>
</div>

<div class="info">
<img src="/images/icons/DATA.png" width="30" />
<p>Data files that are often downloaded from external sources.</p>
</div>

<div class="info">
<img src="/images/icons/BIN.png" width="30" />
<p>An anvi'o bin that describes one or more items. Depending on the context, and item could be a contig, gene cluster, or a user-defined concept.</p>
</div>

<div class="info">
<img src="/images/icons/COLLECTION.png" width="30" />
<p>An anvi'o collection. A collection describes one or more bins.</p>
</div>

<div class="info">
<img src="/images/icons/HMM.png" width="30" />
<p>A collection of HMMs. Besides default anvi'o collections, these can be provided by the user.</p>
</div>

<div class="info">
<img src="/images/icons/JSON.png" width="30" />
<p>A JSON-formatted file.</p>
</div>

<div class="info">
<img src="/images/icons/NEWICK.png" width="30" />
<p>A NEWICK-formatted file.</p>
</div>

<div class="info">
<img src="/images/icons/STATS.png" width="30" />
<p>Summary of numerical information.</p>
</div>

<div class="info">
<img src="/images/icons/SVG.png" width="30" />
<p>An SVG file to generate publication-quality figures.</p>
</div>

<div class="info">
<img src="/images/icons/SUMMARY.png" width="30" />
<p>Anvi'o static summary output HTML page. These are extensive reports of data that can be browsed without an anvi'o installation.</p>
</div>

</div>

<div id="svg"></div>
{% capture network_path %}{{ "/software/anvio/network/network.json" }}{% endcapture %}
{% include _project-anvio-graph.html %}
