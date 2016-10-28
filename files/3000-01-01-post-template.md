---
layout: post
title: "TITLE"
excerpt: "EXCERPT"
modified: 3000-01-01
tags: []
categories: [anvio]
comments: true
authors: []
---

{% capture images %}{{site.url}}/images/anvio/3000-01-01-post-dir{% endcapture %}

{% include _toc.html %}

{:.notice}
{% include _fixthispage.html source="anvio/3000-01-01-template.md" %}

[Link to a post]({% post_url anvio/3000-01-01-template %})

[![Awesome Image]({{images}}/image.png)]({{images}}/image.png){:.center-img .width-80}

{% highlight bash %}
{% endhighlight %}