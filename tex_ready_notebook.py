from html.parser import HTMLParser
import json
import os

os.chdir('./AutomatiqueAvancee')


class MyHTMLParser(HTMLParser):

    attrs_dict = {}

    def handle_starttag(self, tag, attrs):
        for attr in attrs:
            self.attrs_dict[attr[0]] =  attr[1]

def main():
    nb_files = []
    for elem in os.listdir():
        if elem[-6:] == '.ipynb':
            nb_files.append(elem)

    parser = MyHTMLParser()

    for nb_file in nb_files:
        temp_file = nb_file[:-6] + '_temp.ipynb'
        with open(nb_file, 'r') as isf, open(temp_file, 'w') as osf:
            notebook = json.load(isf)
            for cell in notebook['cells']:
                if 'img' in cell.get('source')[0]:
                    parser.feed(cell['source'][0])
                    parser.close()
                    attrs = parser.attrs_dict
                    src = attrs['src']
                    alt = attrs['alt']
                    width = attrs['width']
                    cell['source'] = [
                        r'\begin{figure}[!hbt]', '\n',
                        '\t\t', r'\centering', '\n',
                        '\t\t', r'\includegraphics{{{}}}'.format(src), '\n',
                        '\t\t', r'\caption{{{}}}'.format(alt), '\n',
                        '\t', r'\end{figure}'
                    ]

            json.dump(notebook, osf, indent=2)

if __name__ == "__main__":
    main()
