---
layout: page
title: anvi-script-as-markdown [program]
categories: [anvio]
comments: false
redirect_from: /m/anvi-script-as-markdown
image:
  featurerelative: ../../../images/header.png
  display: true
---

Markdownizides TAB-delmited data with headers in terminal: `cat table.txt | anvi-script-as-markdown`.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/meren.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">A. Murat Eren (Meren)</span><div class="page-author-social-box"><a href="http://meren.org" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:a.murat.eren@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/merenbey" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/meren" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


This program seems to know what its doing. It needs no input material from its user. Good program.


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[markdown-txt](../../artifacts/markdown-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>


## Usage


A helper script to format TAB-delimited files in markdown for bloggers, tutorial writers, those who wish to share example anvi'o outputs as text on GitHub issues, and so on.

Anvi'o programs often generate TAB-delimited files. While this simple format is useful to pass around to other software or share with others, it is not easily interpretable in visual media. The purpose of this script is to make the sharing part simpler for platforms that can render markdown.

You can pipe any TAB-delimited content to this script:

```
cat file.txt | anvi-sript-as-markdown
```

## Examples

Assume a TAB-delmited file with many lines:

``` bash
wc -l additional_view_data.txt

301 additional_view_data.txt
```

Contents of which look like this:

``` bash
head -n 10 additional_view_data.txt

contig	categorical_1	categorical_2	text_layer_01	numerical	bars_main!A	bars_main!B	bars_main!C
backrest	b	y	nmwje	2.78	278	23	1
backward	b	x	bqmyujr psrd doefhi	2.49	249	52	2
backwind	b	y	hkfer lchpmzix	2.69	269	32	3
backyard	b	x	advoe bfkyhmg	2.05	205	96	4
bacteria	b	x	lqmcwn hywco	2.63	263	38	5
bacterin	b		vxqdmn	2.98	298	3	6
baetylus	b	x	fkgpydi owgyhfx xwlpj	2.19	219	82	7
bagpiped	b	y	ijmnur	2.12	212	89	8
balconet	b	y	ecizgs	2.89	289	12	9
```

### Default run

<div class="codeblock" markdown="1">
head &#45;n 10 additional_view_data.txt | <span class="artifact&#45;n">[anvi&#45;script&#45;as&#45;markdown](/software/anvio/help/main/programs/anvi&#45;script&#45;as&#45;markdown)</span>
</div>

which is rendered as,

|**contig**|**categorical_1**|**categorical_2**|**text_layer_01**|**numerical**|**bars_main!A**|**bars_main!B**|**bars_main!C**|
|:--|:--|:--|:--|:--|:--|:--|:--|
|backrest|b|y|nmwje|2.78|278|23|1|
|backward|b|x|bqmyujr psrd doefhi|2.49|249|52|2|
|backwind|b|y|hkfer lchpmzix|2.69|269|32|3|
|backyard|b|x|advoe bfkyhmg|2.05|205|96|4|
|bacteria|b|x|lqmcwn hywco|2.63|263|38|5|
|bacterin|b||vxqdmn|2.98|298|3|6|
|baetylus|b|x|fkgpydi owgyhfx xwlpj|2.19|219|82|7|
|bagpiped|b|y|ijmnur|2.12|212|89|8|
|balconet|b|y|ecizgs|2.89|289|12|9|

### Limit the number of lines shown

``` bash
cat additional_view_data.txt | anvi-script-as-markdown --max-num-lines-to-show 10
```

which is rendered as,

|**contig**|**categorical_1**|**categorical_2**|**text_layer_01**|**numerical**|**bars_main!A**|**bars_main!B**|**bars_main!C**|
|:--|:--|:--|:--|:--|:--|:--|:--|
|backrest|b|y|nmwje|2.78|278|23|1|
|backward|b|x|bqmyujr psrd doefhi|2.49|249|52|2|
|backwind|b|y|hkfer lchpmzix|2.69|269|32|3|
|backyard|b|x|advoe bfkyhmg|2.05|205|96|4|
|bacteria|b|x|lqmcwn hywco|2.63|263|38|5|
|bacterin|b||vxqdmn|2.98|298|3|6|
|baetylus|b|x|fkgpydi owgyhfx xwlpj|2.19|219|82|7|
|bagpiped|b|y|ijmnur|2.12|212|89|8|
|balconet|b|y|ecizgs|2.89|289|12|9|
|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|

### Code columns

``` bash
cat additional_view_data.txt | anvi-script-as-markdown --max-num-lines-to-show 10 \
                                                       --code-column contig
```

which is rendered as,

|**contig**|**categorical_1**|**categorical_2**|**text_layer_01**|**numerical**|**bars_main!A**|**bars_main!B**|**bars_main!C**|
|:--|:--|:--|:--|:--|:--|:--|:--|
|`backrest`|b|y|nmwje|2.78|278|23|1|
|`backward`|b|x|bqmyujr psrd doefhi|2.49|249|52|2|
|`backwind`|b|y|hkfer lchpmzix|2.69|269|32|3|
|`backyard`|b|x|advoe bfkyhmg|2.05|205|96|4|
|`bacteria`|b|x|lqmcwn hywco|2.63|263|38|5|
|`bacterin`|b||vxqdmn|2.98|298|3|6|
|`baetylus`|b|x|fkgpydi owgyhfx xwlpj|2.19|219|82|7|
|`bagpiped`|b|y|ijmnur|2.12|212|89|8|
|`balconet`|b|y|ecizgs|2.89|289|12|9|
|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|

### Exclude columns from the output

``` bash
cat additional_view_data.txt | anvi-script-as-markdown --max-num-lines-to-show 10 \
                                                       --code-column contig \
                                                       --exclude-columns 'bars_main!A,bars_main!B,bars_main!C'
```

which is rendered as,

|**contig**|**categorical_1**|**categorical_2**|**text_layer_01**|**numerical**|
|:--|:--|:--|:--|:--|
|`backrest`|b|y|nmwje|2.78|
|`backward`|b|x|bqmyujr psrd doefhi|2.49|
|`backwind`|b|y|hkfer lchpmzix|2.69|
|`backyard`|b|x|advoe bfkyhmg|2.05|
|`bacteria`|b|x|lqmcwn hywco|2.63|
|`bacterin`|b||vxqdmn|2.98|
|`baetylus`|b|x|fkgpydi owgyhfx xwlpj|2.19|
|`bagpiped`|b|y|ijmnur|2.12|
|`balconet`|b|y|ecizgs|2.89|
|(...)|(...)|(...)|(...)|(...)|


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-script-as-markdown.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-script-as-markdown) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
