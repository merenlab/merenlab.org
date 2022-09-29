# -*- coding: utf-8 -*-
# an ugly hack to convert some stuff into other stuff...


# EDIT THESE #####################################################################
names_to_highlight = {'Eren AM': None,
                      'Delmont TO': range(2015, 2020),
                      'Esen ÖC': range(2015, 2021),
                      'Esen OC': range(2015, 2021),
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

import sys
from datetime import datetime

try:
    import anvio.utils as u
    from anvio.errors import ConfigError
except:
    sys.stderr.write("This program requires anvi'o to be installed :/\n")
    sys.exit(-1)


class Publications:
    def __init__(self, pubs_file_path='_data/pubs.yaml'):
        """Takes a YAML file, generates pubs output"""

        self.info = {}

        self.pubs_dict = {}
        self.journals_list = []
        self.authors_list = []
        self.recent_authors_list = []
        self.author_links = {}

        self.pubs_file_path = pubs_file_path


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


    def get_abbreviated_name_from_full_name(self, author_name):
        """Takes full name (First Middle Last), return abbreviated name (Last FM)"""

        names = author_name.replace('.', '').split()

        abbreviated_name = f"{names[-1]} {''.join([n[0] for n in names[:-1]])}"

        if '*' in abbreviated_name:
            abbreviated_name = abbreviated_name.replace('*', '') + '*'

        if '+' in abbreviated_name:
            abbreviated_name = abbreviated_name.replace('+', '') + '+'

        return abbreviated_name


    def parse_pubs_txt(self):
        self.pubs = u.get_yaml_as_dict(self.pubs_file_path)

        # {'doi': '10.1101/2020.11.01.361691',
        #  'title': 'The genetic and ecological landscape of plasmids in the human gut',
        #  'authors': 'Michael K Yu*, Emily C Fogarty*, A. Murat Eren',
        #  'journal': 'bioRxiv',
        #  'year': 2022,
        #  'volume': None,
        #  'number': None,
        #  'pages': None,
        #  'additional_info': {
        #       'highlights': ['A study that aims to shed light on <b>the ecology and evolution of one of the most critical yet poorly studied aspects of microbial life -- naturally occurring plasmids</b>.',
        #                      'Uses state-of-the-art machine learning strategies to identify <b>over 60,000 plasmids</b> from human gut metagenomes, which represents a <b>200-fold increase</b> in the number of known plasmids to date that were detectable in healthy humans.',
        #                      "Defines hundreds of '<b>plasmid systems</b>', and demonstrates that naturally occurring plasmids are not static entities, but <b>their evolution is driven by the need to respond to the environment, and their ecology cannot be simply explained by bacterial taxonomy and distribution patterns of their putative hosts</b>."],
        #       'featured_image': '/images/pubs/plasmid_systems.png'
        #  }
        # }

        for pub in self.pubs:
            pub['co_first_authors'] = []
            pub['co_senior_authors'] = []

            # turn author names from "FIRST M LAST" form to "LAST FM" form.
            pub['authors'] = [self.get_abbreviated_name_from_full_name(a.strip()) for a in pub['authors'].split(',')]

            for author_name in pub['authors']:
                if author_name.endswith('*'):
                    pub['co_first_authors'].append(author_name[:-1])
                elif  author_name.endswith('+'):
                    pub['co_senior_authors'].append(author_name[:-1])

            pub['authors'] = [a[:-1] if (a.endswith('*') or a.endswith('+')) else a for a in pub['authors']]

            if pub['volume'] and pub['number'] and pub['pages']:
                pub['issue'] = f"{pub['volume']}({pub['number']}):{pub['pages']}"
            elif pub['volume'] and pub['number']:
                pub['issue'] = f"{pub['volume']}({pub['number']})"
            elif pub['volume'] and pub['pages']:
                pub['issue'] = f"{pub['volume']}:{pub['pages']}"
            elif pub['volume']:
                pub['issue'] = f"{pub['volume']}"
            else:
                pub['issue'] = None

            year = pub['year']
            if year not in self.pubs_dict:
                self.pubs_dict[year] = [pub]
            else:
                self.pubs_dict[year].append(pub)


    def get_markdown_text_for_pub(self, pub):
        pub_md = []

        A = lambda s: pub_md.append(s)

        if 'read_link' in pub:
            read_link = pub['read_link']
        else:
            read_link = f"https://doi.org/{pub['doi']}"

        A(f'''<a id="{pub['doi']}">&nbsp;</a>''')
        A('<div class="pub">')
        A(f'''<div class='altmetric-embed' data-badge-type='donut' data-doi="{pub['doi']}"></div>''')
        A(f'''<div class="__dimensions_badge_embed__" data-doi="{pub['doi']}" data-hide-zero-citations="true" data-legend="hover-bottom" data-style="small_circle"></div>''')
        A(f'''    <span class="pub-title"><a href="{read_link}" target="_new">{pub['title']}</a></span>''')
        A(f'''    <span class="pub-authors">{self.get_author_highlights(pub, pub['year'])}</span>''')

        # take care of co-first / co-senior authors
        if pub['co_first_authors'] and not pub['co_senior_authors']:
            A('    <span class="pub-co-first-authors"><sup>☯</sup>Co-first authors</span>')
        elif pub['co_first_authors'] and pub['co_senior_authors']:
            A('    <span class="pub-co-first-authors"><sup>☯</sup>Co-first authors; <sup>‡</sup>Co-senior authors</span>')
        elif pub['co_senior_authors'] and not pub['co_first_authors']:
            A('    <span class="pub-co-first-authors"><sup>‡</sup>Co-senior authors</span>')

        # add the publication highlights:
        if 'additional_info' in pub and pub['additional_info']['highlights']:
            I = pub['additional_info']
            A('    <div class="%s">' % ('pub-info' if pub['additional_info']['featured_image'] else 'pub-info-no-image'))

            if I['featured_image']:
                A('    <div class="pub-featured-image">')
                A('    <a href="%s"><img src="%s" style="max-width: 100px; max-height: 80px; width: auto; border: none; height: auto; margin: 0 auto; display: block; transform: translateY(15%%);"/></a>' % (I['featured_image'], I['featured_image']))
                A('    </div>')

            if I['highlights']:
                A('    <div class="%s">' % ('pub-highlights' if I['featured_image'] else 'pub-highlights-no-image'))
                A('    %s' % '<br>'.join(['<span style="display: inline-block; padding-bottom: 5px;">- %s</span>' % h for h in I['highlights']]))
                A('    </div>')

            A('    </div>')


        scholar_link = f'''http://scholar.google.com/scholar?hl=en&q={pub['title'].replace(' ', '+')}'''
        additional_links = f'''| 🔍 <a href="{scholar_link}" target="_blank">Google Scholar</a> | 🔗 <a href="https://doi.org/{pub['doi']}" target="_blank">doi:{pub['doi']}</a>'''

        if pub['issue']:
            A(f'''    <span class="pub-journal"> 📚 <b>{pub['journal']}</b>, {pub['issue']} {additional_links}</span>''')
        else:
            A(f'''    <span class="pub-journal"> 📚 <b>{pub['journal']}</b> {additional_links}</span>''')
        A('</div>\n')

        return '\n'.join(pub_md)


    def store_markdown_output_for_pubs(self, output_file_path):
        #years = ''.join(['<a href="#%s"><span class="category-item">%s</span></a>' % (y, y) for y in sorted(list(self.pubs_dict.keys()), reverse=True)])

        output_file = open(output_file_path, 'w')
        W = lambda s: output_file.write(s + '\n')

        W('---')
        W('layout: publications')
        W('modified: %s' % datetime.today().strftime('%Y-%m-%d'))
        W('comments: false')
        W('image:')
        W('   display: true')
        W('   feature: eel-pond.jpg')
        W('---\n')

        W('''<script type='text/javascript' src='https://d1bxh8uas1mnw7.cloudfront.net/assets/embed.js'></script>\n''')
        W('''<script async src="https://badge.dimensions.ai/badge.js" charset="utf-8"></script>\n''')

        #W('<div class="category-box">\n%s\n</div>\n' % years)

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
