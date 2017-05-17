---
layout: post
title: "The new anvi'o 'news' tab"
excerpt: "Is the big brother watching you? If he does, what does he see?"
modified: 2017-05-16
tags: [privacy, interface]
categories: [anvio]
comments: true
authors: [meren]
---


{% include _toc.html %}

If you are using anvi'o `v2.3.2` or later, the anvi'o interactive interface will have a 'news' tab for you on the right-hand side of your screen.

The purpose of this post is to clarify some points regarding the news tab, and create a resource under which we can address your further questions or concerns.

## Why a news tab?

It has always been very hard for us to communicate with our users. For instance, we have been keeping our users up-to-date through our blog on our lab web page, anvi'o discussion group, or through our Twitter accounts. However, we realize that there is a lot going on in all these places, and it is easy to miss things that may be important to an anvi'o user. This could be a new release, a workshop, or a critical bug that requires an update. The purpose of our new tab is to address this need.

## Should I be concerned about my privacy?

No. Anvi'o does not collect any kind of data from you. The anvi'o interactive interface you are using with `v2.3.2` or later every once in awhile your browser accesses to [this file](https://github.com/merenlab/anvio/blob/master/NEWS.md) on GitHub through a simple [HTTP](https://en.wikipedia.org/wiki/Hypertext_Transfer_Protocol) request. If there is anything new in that file, it lights up the 'News' tab. To do this, it stores two [HTTP cookie](https://en.wikipedia.org/wiki/HTTP_cookie)s in your browser. One of these cookies keep track of when does the news was last checked, so you send a request only once in a day. The second cookie keeps the signature (hash sum) of the news file at the time of access, so if nothing has chanced the interface does not SPAM you with fake news alerts.

Because you are sending a read request to a file that is publicly accessible on GitHub, your IP address or location is not accessible to anyone, including the anvi'o developers.

## How often?

We will use this feature only to keep our users up-to-date regarding to anvi'o related developments. 

---

Please let us know about your questions, suggestions, or concerns by leaving a comment or sending an e-mail. Thank you!

<div style="margin:100px">&nbsp;</div>