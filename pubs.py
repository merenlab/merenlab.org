# -*- coding: utf-8 -*-
# an ugly hack to convert some stuff into other stuff...

# people who have links
author_links = {}

bpc_people = ['Eren, A. M', 'Delmont, T. O.', 'Esen, Ö. C.']

keep_pubs_after = 2009

pubs_dict = {}
journals_list = []
authors_list = []
recent_authors_list = []

# takes an EndNote library exported as a TXT file. here is a sample line from this txt file:
#
# Winterberg, K. M., and Reznikoff, W. S. (2007). "Screening transposon mutant libraries using full-genome oligonucleotide microarrays." Methods Enzymol, 421, 110-25.
#
f = open('pubs.txt')

def get_author_links(authors_str):
    for author in bpc_people:
        authors_str = authors_str.replace(author, '<span class="pub-member-author">%s</span>' % (author))

    return authors_str


for line in [l.strip() for l in f.readlines()]:
    if line.find('(ed.)') > 0 or line.find('(eds.)') > 0:
        continue
    p_s = line.find(' (')
    p_e = p_s + 6
    if not p_s > 0:
        continue
    if not line[p_e] == ')':
        continue

    year = int(line[p_s + 2:p_e])

    if year < keep_pubs_after:
        continue

    authors = line[0: p_s]

    q_s = line.find(' "', p_e)
    if not q_s > 0:
        continue
    q_e = line.find('."', q_s)

    if not q_e > 0:
        q_e = line.find('?"', q_s)
        if not q_e > 0:
            continue

    title = line[q_s + 2:q_e + 1]

    c = line.find(', ', q_e + 2)
    if not c > 0:
        continue

    journal = line[q_e + 3:c]

    issue = line[c + 2:-1]


    # ad hoc fixes for journal names
    journal = journal.replace('The ISME journal', 'ISME J')
    journal = journal.replace('Proceedings of the National Academy of Sciences of the United States of America', 'Proc Natl Acad Sci U S A')
    journal = journal.replace('Proceedings of the National Academy of Sciences', 'Proc Natl Acad Sci U S A')
    journal = journal.replace('Frontiers in Microbiology', 'Front Microbiol')
    journals_list.append(journal)

    authors = authors.replace('Esen, Ö.,', 'Esen, Ö. C.,')
    authors = authors.replace('Murat Eren, A.,', 'Eren, A. M.,')

    if not pubs_dict.has_key(year):
        pubs_dict[year] = [{'authors': authors, 'title': title, 'journal': journal, 'issue': issue}]
    else:
        pubs_dict[year].append({'authors': authors, 'title': title, 'journal': journal, 'issue': issue})


    if authors.count(',') == 1:
        authors_list.append(authors)
        if year > 2004:
            recent_authors_list.append(authors)
    else:
        for author in [a + '.' if not a.endswith('.') else a for a in authors.replace('and ', '').split('., ')]:
            authors_list.append(author)
            if year > 2004:
                recent_authors_list.append(author)


years = ''.join(['<a href="#%s"><span class="category-item">%s <small>(%d)</small></span></a>' % (y, y, len(pubs_dict[y])) for y in sorted(pubs_dict.keys(), reverse=True)])

top_journals = ", ".join(['<b>%s</b> (<i>%d</i>)' % (x[1], x[0]) for x in sorted([(journals_list.count(journal), journal) for journal in set(journals_list)], reverse = True)[0:25]])

print """---
layout: publications
modified: 2015-02-05
comments: false
---

"""
# print "<h1>Journals</h1>"
# print top_journals
# print
print '<div class="category-box">'
print years
print '</div>'
print

for year in sorted(pubs_dict.keys(), reverse=True):
    print '<a name="%s">&nbsp;</a>' % year
    print '<h1>%s</h1>' % year
    print
    for pub in pubs_dict[year]:
        print '<div class="pub">'
        print '    <h3><a href="http://scholar.google.com/scholar?hl=en&q=%s" target="_new">%s</a></h3>' % ('http://scholar.google.com/scholar?hl=en&q=%s' % (pub['title'].replace(' ', '+')), pub['title'])
        print '    <span class="pub-authors">%s</span>' % get_author_links(pub['authors'])
        print '    <span class="pub-journal"><i>%s</i>. <b>%s</b></span>' % (pub['journal'], pub['issue'])
        print '</div>'
        print
    print

