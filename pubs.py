# -*- coding: utf-8 -*-
# an ugly hack to convert some stuff into other stuff...


# EDIT THESE #####################################################################
names_to_highlight = {'Eren AM': None,
                      'Delmont TO': range(2015, 2020),
                      'Esen ÖC': range(2015, 2021),
                      'Yu MK': None,
                      'Lee STM': None,
                      'Shaiber A': None,
                      'Kiefl E': None,
                      'Cui S': None,
                      'Watson AR': None,
                      'Lolans K': None,
                      'Schmid AC': None,
                      'Yousef M': None,
                      'Veseli I': None,
                      'Miller SE': None,
                      'Schechter MS': None,
                      'Fink I': None,
                      'Pan JN': None,
                      'Yousef M': None,
                      'Fogarty EC': None,
                      'Trigodet F': None,}

journal_name_fixes = [('The ISME journal', 'ISME J'),
                      ('Proceedings of the National Academy of Sciences of the United States of America', 'Proc Natl Acad Sci U S A'),
                      ('Proceedings of the National Academy of Sciences', 'Proc Natl Acad Sci U S A'),
                      ('Frontiers in Microbiology', 'Front Microbiol')]

keep_pubs_after_year = 2009
##################################################################################

import os
import sys
from datetime import datetime

try:
    import anvio.utils as u
    from anvio.errors import ConfigError
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


    def get_author_highlights(self, pub, year):
        authors_str = []
        for author in pub['authors']:
            if author in pub['co_first_authors']:
                author_h = author + '<sup>☯</sup>'
            elif author in pub['co_senior_authors']:
                author_h = author + '<sup>‡</sup>'
            else:
                author_h = author

            if author in names_to_highlight:
                if not names_to_highlight[author]:
                    authors_str.append('<span class="pub-member-author">%s</span>' % (author_h))
                elif int(year) in names_to_highlight[author]:
                    authors_str.append('<span class="pub-member-author">%s</span>' % (author_h))
                else:
                    authors_str.append(author_h)
            else:
                authors_str.append(author_h)

        return ', '.join(authors_str)


    def parse_pubs_txt(self):
        if os.path.exists(self.pubs_info_file_path):
            self.info = u.get_TAB_delimited_file_as_dictionary(self.pubs_info_file_path)

        pubs_header = u.get_columns_of_TAB_delim_file(self.pubs_file_path, include_first_column=True)
        headers_expected = ['Authors', 'Title', 'Publication', 'Volume', 'Number', 'Pages', 'Year', 'doi']
        missing_headers = [h for h in pubs_header if h not in headers_expected]
        if len(missing_headers):
            raise ConfigError("Sorry, the pubs.txt seems to be missing some of the headers that are mandatory. Each of \
                               the columns in the following list must be present in this file: %s (hint: yours do not have\
                               the following: %s)." % (', '.join(headers_expected), ', '.join(missing_headers)))

        self.pubs_txt = u.get_TAB_delimited_file_as_dictionary(self.pubs_file_path, indexing_field=pubs_header.index('doi'))

        for doi in self.pubs_txt:
            authors = []
            co_first_authors = []
            co_senior_authors = []
            p = self.pubs_txt[doi]

            for author in [_.strip() for _ in p['Authors'].split(';')]:
                if not len(author):
                    continue

                author_last_name, author_first_name_raw = [_.strip() for _ in author.split(',')]
                author_first_name = ''.join([n[0] for n in author_first_name_raw.split()])
                author_final_name = '%s %s' % (author_last_name, author_first_name)

                if author_first_name_raw.endswith('*'):
                    co_first_authors.append(author_final_name)
                elif  author_first_name_raw.endswith('+'):
                    co_senior_authors.append(author_final_name)

                authors.append(author_final_name)

            if p['Number']:
                issue = '%s(%s):%s' % (p['Volume'], p['Number'], p['Pages'])
            elif p['Volume'] and p['Pages']:
                issue = '%s:%s' % (p['Volume'], p['Pages'])
            else:
                issue = None

            year = p['Year'].strip()
            pub_entry = {'authors': authors, 'title': p['Title'], 'journal': p['Publication'], 'issue': issue, 'doi': doi, 'year': year, 'co_first_authors': co_first_authors, 'co_senior_authors': co_senior_authors}

            if year not in self.pubs_dict:
                self.pubs_dict[year] = [pub_entry]
            else:
                self.pubs_dict[year].append(pub_entry)


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
        A('''<div class="__dimensions_badge_embed__" data-doi="%s" data-hide-zero-citations="true" data-legend="hover-bottom" data-style="small_circle"></div>''' % pub['doi'])
        if pub['doi']:
            A('    <span class="pub-title"><a href="%s" target="_new">%s</a></span>' % (' https://doi.org/%s' % (pub['doi']), pub['title']))
        else:
            A('    <span class="pub-title"><a href="http://scholar.google.com/scholar?hl=en&q=%s" target="_new">%s</a></span>' % ('http://scholar.google.com/scholar?hl=en&q=%s' % (pub['title'].replace(' ', '+')), pub['title']))
        A('    <span class="pub-authors">%s</span>' % self.get_author_highlights(pub, pub['year']))

        if pub['co_first_authors'] and not pub['co_senior_authors']:
            A('    <span class="pub-co-first-authors"><sup>☯</sup>Co-first authors</span>')
        elif pub['co_first_authors'] and pub['co_senior_authors']:
            A('    <span class="pub-co-first-authors"><sup>☯</sup>Co-first authors; <sup>‡</sup>Co-senior authors</span>')
        elif pub['co_senior_authors'] and not pub['co_first_authors']:
            A('    <span class="pub-co-first-authors"><sup>‡</sup>Co-senior authors</span>')

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

        if pub['issue']:
            A('    <span class="pub-journal"><b>%s</b>, %s <a href="https://doi.org/%s" target="_blank">🔗</a></span>' % (pub['journal'], pub['issue'], pub['doi']))
        else:
            A('    <span class="pub-journal"><b>%s</b> <a href="https://doi.org/%s" target="_blank">🔗</a></span>' % (pub['journal'], pub['doi']))
        A('</div>\n')

        return '\n'.join(pub_md)


    def store_markdown_output_for_pubs(self, output_file_path):
        # years = ''.join(['<a href="#%s"><span class="category-item">%s <small>(%d)</small></span></a>' % (y, y, len(self.pubs_dict[y])) for y in sorted(list(self.pubs_dict.keys()), reverse=True)])
        years = ''.join(['<a href="#%s"><span class="category-item">%s</span></a>' % (y, y) for y in sorted(list(self.pubs_dict.keys()), reverse=True)])

        output_file = open(output_file_path, 'w')
        W = lambda s: output_file.write(s + '\n')

        W('---')
        W('layout: publications')
        W('modified: %s' % datetime.today().strftime('%Y-%m-%d'))
        W('comments: false')
        W('---\n')

        W('''<script type='text/javascript' src='https://d1bxh8uas1mnw7.cloudfront.net/assets/embed.js'></script>\n''')
        W('''<script async src="https://badge.dimensions.ai/badge.js" charset="utf-8"></script>\n''')

        W('<div class="category-box">\n%s\n</div>\n' % years)

        W('{:.notice}\n')
        W("This page lists publications that are most reflective of our interests. For a complete list, please see <a href='https://scholar.google.com/citations?user=GtLLuxoAAAAJ&view_op=list_works&sortby=pubdate' target='_blank'>Meren's Google Scholar page</a>.\n")

        for year in sorted(list(self.pubs_dict.keys()), reverse=True):
            W('## %s\n' % (year))

            for pub in self.pubs_dict[year]:
                W(self.get_markdown_text_for_pub(pub))

            W('')


if __name__ == '__main__':
    pubs = Publications()
    try:
        pubs.parse_pubs_txt()
        pubs.store_markdown_output_for_pubs('publications/index.md')
    except ConfigError as e:
        print(e)
        sys.exit(-1)
