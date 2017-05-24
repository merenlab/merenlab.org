# -*- coding: utf-8 -*-
# an ugly hack to convert some stuff into other stuff...


# EDIT THESE #####################################################################
names_to_highlight = ['Eren, A. M', 
                      'Delmont, T. O.',
                      'Esen, Ö. C.',
                      'Lee, S. T. M.',
                      'Shaiber, A.',
                      'Kiefl, E.']

journal_name_fixes = [('The ISME journal', 'ISME J'),
                      ('Proceedings of the National Academy of Sciences of the United States of America', 'Proc Natl Acad Sci U S A'),
                      ('Proceedings of the National Academy of Sciences', 'Proc Natl Acad Sci U S A'),
                      ('Frontiers in Microbiology', 'Front Microbiol')]

keep_pubs_after_year = 2009
##################################################################################

import os
import sys

try:
    import anvio.utils as u
except:
    sys.stderr.write("This program requires anvi'o to be installed :/\n")
    sys.exit(-1)


class Publications:
    def __init__(self, pubs_file_path='pubs.txt', pubs_info_file_path='pubs_info.txt'):
        """Takes an EndNote library exported a TXT file (`pubs_file_path`), and an optional\
           TAB-delimited info file path with DOI identifiers (`pubs_info_file_path`), and\
           generates some Markdown formatted output.

           Here is an info line from the EndNote:

                Winterberg, K. M., and Reznikoff, W. S. (2007). "Screening transposon mutant libraries using full-genome oligonucleotide microarrays." Methods Enzymol, 421, 110-25.

           Absolute matching to this format is required.

           Expected headers in the TAB-delimited pubs info file are 'doi', 'highlights',\
           and 'featured_image'.

                - doi: The DOI of the pub matching to a pubs file path entry.
                - highlights: Brief bullet points about the work. Each pont must be separated\
                              from the rest with a ';' character. HTML tags are OK.
                - featured_image: A URL to an image.

           If things are not working, feel free to write to meren at uchicago.edu
        """

        self.info = {}

        self.pubs_dict = {}
        self.journals_list = []
        self.authors_list = []
        self.recent_authors_list = []
        self.author_links = {}

        self.pubs_file_path = pubs_file_path
        self.pubs_info_file_path = pubs_info_file_path


    def get_author_highlights(self, authors_str):
        for author in names_to_highlight:
            authors_str = authors_str.replace(author, '<span class="pub-member-author">%s</span>' % (author))

        return authors_str


    def parse_bookends_input(self):
        bad_entries = []

        if os.path.exists(self.pubs_info_file_path):
            self.info = u.get_TAB_delimited_file_as_dictionary(self.pubs_info_file_path)

        for line in [l.strip() for l in open(self.pubs_file_path, 'rU').readlines()]:
            if line.find('(ed.)') > 0 or line.find('(eds.)') > 0:
                bad_entries.append((line, 'ed/eds. found...'), )
                continue

            p_s = line.find(' (')
            p_e = p_s + 6
            if not p_s > 0:
                bad_entries.append((line, 'p_s <= 0...'), )
                continue
            if not line[p_e] == ')':
                bad_entries.append((line, 'p_e != )...'), )
                continue

            doi = None
            if line.split()[-1].strip().startswith('doi:'):
                doi = line.split('doi:')[1].strip()
                line = line.split('doi:')[0].strip()

            year = int(line[p_s + 2:p_e])

            if year < keep_pubs_after_year:
                bad_entries.append((line, 'year >= keep_pubs_after...'), )
                continue

            authors = line[0: p_s]

            q_s = line.find(' "', p_e)
            if not q_s > 0:
                bad_entries.append((line, 'q_s <= 0...'), )
                continue
            q_e = line.find('."', q_s)

            if not q_e > 0:
                q_e = line.find('?"', q_s)
                if not q_e > 0:
                    bad_entries.append((line, 'q_e <= 0...'), )
                    continue

            title = line[q_s + 2:q_e + 1]

            c = line.find(', ', q_e + 2)
            if not c > 0:
                bad_entries.append((line, 'c <= 0...'), )
                continue

            journal = line[q_e + 3:c]

            issue = line[c + 2:-1]

            # ad hoc fixes for journal names
            for bad_form, good_form in journal_name_fixes:
                journal = journal.replace(bad_form, good_form)

            self.journals_list.append(journal)

            authors = authors.replace('Esen, Ö.,', 'Esen, Ö. C.,')
            authors = authors.replace('Murat Eren, A.,', 'Eren, A. M.,')

            if year not in self.pubs_dict:
                self.pubs_dict[year] = [{'authors': authors, 'title': title, 'journal': journal, 'issue': issue, 'doi': doi, 'year': year}]
            else:
                self.pubs_dict[year].append({'authors': authors, 'title': title, 'journal': journal, 'issue': issue, 'doi': doi, 'year': year})

            if authors.count(',') == 1:
                self.authors_list.append(authors)
                if year > 2004:
                    self.recent_authors_list.append(authors)
            else:
                for author in [a + '.' if not a.endswith('.') else a for a in authors.replace('and ', '').split('., ')]:
                    self.authors_list.append(author)
                    if year > 2004:
                        self.recent_authors_list.append(author)


        # check for failed entries
        if len(bad_entries):
            print("Some entries failed. Quitting.")
            print()
            for tpl in bad_entries:
                print(' - Failed (reason: "%s"): %s' % (tpl[1], tpl[0]))

            sys.exit(-2)


    def get_markdown_text_for_pub(self, pub):
        """Gets a dictionary `pub`, returns a markdown formatted text.

           An example pub:

                {'authors': 'McLellan, S. L., and Eren, A. M.',
                 'doi': '10.1016/j.tim.2014.08.002',
                 'issue': '22(12), 697-706',
                 'title': 'Discovering new indicators of fecal pollution.',
                 'journal': 'Trends Microbiol',
                 'year': 2014}
        """

        pub_md = []

        A = lambda s: pub_md.append(s)

        A('<div class="pub">')
        A('''<div class='altmetric-embed' data-badge-type='donut' data-doi="%s"></div>''' % pub['doi'])
        if pub['doi']:
            A('    <h3><a href="%s" target="_new">%s</a></h3>' % (' https://doi.org/%s' % (pub['doi']), pub['title']))
        else:
            A('    <h3><a href="http://scholar.google.com/scholar?hl=en&q=%s" target="_new">%s</a></h3>' % ('http://scholar.google.com/scholar?hl=en&q=%s' % (pub['title'].replace(' ', '+')), pub['title']))
        A('    <span class="pub-authors">%s</span>' % self.get_author_highlights(pub['authors']))

        if pub['doi'] in self.info:
            info = self.info[pub['doi']]
            A('    <div class="%s">' % ('pub-info' if info['featured_image'] else 'pub-info-no-image'))

            if info['featured_image']:
                A('    <div class="pub-featured-image">')
                A('    <a href="%s"><img src="%s" style="max-width: 100px; max-height: 80px; width: auto; border: none; height: auto; margin: 0 auto; display: block; transform: translateY(15%%);"/></a>' % (info['featured_image'], info['featured_image']))
                A('    </div>')

            highlights = info['highlights'].split(';') if info['highlights'] else None
            if highlights:
                A('    <div class="%s">' % ('pub-highlights' if info['featured_image'] else 'pub-highlights-no-image'))
                A('    %s' % '<br>'.join(['<span style="display: inline-block; padding-bottom: 5px;">- %s</span>' % h for h in highlights]))
                A('    </div>')

            A('    </div>')

        A('    <span class="pub-journal"><i>%s</i>. <b>%s</b></span>' % (pub['journal'], pub['issue']))
        A('</div>\n')

        return '\n'.join(pub_md)


    def store_markdown_output_for_pubs(self, output_file_path, include_top_journals = False):
        years = ''.join(['<a href="#%s"><span class="category-item">%s <small>(%d)</small></span></a>' % (y, y, len(self.pubs_dict[y])) for y in sorted(list(self.pubs_dict.keys()), reverse=True)])
        top_journals = ", ".join(['<b>%s</b> (<i>%d</i>)' % (x[1], x[0]) for x in sorted([(self.journals_list.count(journal), journal) for journal in set(self.journals_list)], reverse = True)[0:25]])

        output_file = open(output_file_path, 'w')
        W = lambda s: output_file.write(s + '\n')

        W('---')
        W('layout: publications')
        W('modified: 2015-02-05')
        W('comments: false')
        W('---\n')

        W('''<script type='text/javascript' src='https://d1bxh8uas1mnw7.cloudfront.net/assets/embed.js'></script>\n''')

        if include_top_journals:
            W("<h1>Journals</h1>\n%s\n" % top_journals)

        W('<div class="category-box">\n%s\n</div>\n' % years)

        for year in sorted(list(self.pubs_dict.keys()), reverse=True):
            W('<a name="%s">&nbsp;</a>' % year)
            W('<h1>%s</h1>\n' % year)

            for pub in self.pubs_dict[year]:
                W(self.get_markdown_text_for_pub(pub))

            W('')


if __name__ == '__main__':
    pubs = Publications()
    pubs.parse_bookends_input()
    pubs.store_markdown_output_for_pubs('publications/index.md')
