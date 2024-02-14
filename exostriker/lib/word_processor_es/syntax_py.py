# syntax.py
# adapted from = https://wiki.python.org/moin/PyQt/Python%20syntax%20highlighting

import sys

from PyQt6.QtCore import QRegularExpression as QRegExp
from PyQt6.QtGui import QColor, QTextCharFormat, QFont, QSyntaxHighlighter

quote = "%s%s%s" % (chr(39), chr(39), chr(39))
dquote = "%s%s%s" % (chr(34), chr(34), chr(34))

def format(color, style=''):
    '''Return a QTextCharFormat with the given attributes.
    '''
    _color = QColor()
    _color.setNamedColor(color)

    _format = QTextCharFormat()
    _format.setForeground(_color)
    if 'bold' in style:
        _format.setFontWeight(QFont.Weight.Bold)
    if 'italic' in style:
        _format.setFontItalic(True)
    if 'italicbold' in style:
        _format.setFontItalic(True)
        _format.setFontWeight(QFont.Weight.Bold)
    return _format

mybrawn = ("#7E5916")
# Syntax styles that can be shared by all languages
STYLES = {
    'keyword': format('#2C2CC8', 'bold'),
    'operator': format('darkred'),
    'brace': format('darkred'),
    'defclass': format('#cc0000', 'bold'),
    'classes': format('#cc0000', 'bold'),
    'Qtclass': format('black', 'bold'),
    'string': format(mybrawn),
    'string2': format('#42923b', 'italic'),
    'comment': format('#42923b', 'italic'),
    'self': format('#D63030', 'italicbold'),
    'selfnext': format('#2e3436', 'bold'),
    'Qnext': format('#2e3436', 'bold'),
    'numbers': format('#C82C2C'),
}

class Highlighter(QSyntaxHighlighter):
    '''Syntax highlighter for the Python language.
    '''
    # Python keywords
    keywords = [
        'and', 'assert', 'break', 'class', 'continue', 'def',
        'del', 'elif', 'else', 'except', 'exec', 'finally',
        'for', 'from', 'global', 'if', 'import', 'in',
        'is', 'lambda', 'not', 'or', 'pass', 'print',
        'raise', 'return', 'super', 'try', 'while', 'yield',
        'None', 'True', 'False',
    ]

    # Python operators
    operators = [
        '=',
        # Comparison
        '==', '!=', '<', '<=', '>', '>=',
        # Arithmetic
        '\+', '-', '\*', '/', '//', '\%', '\*\*',
        # In-place
        '\+=', '-=', '\*=', '/=', '\%=',
        # Bitwise
        '\^', '\|', '\&', '\~', '>>', '<<',
    ]

    # Python braces
    braces = [
        '\{', '\}', '\(', '\)', '\[', '\]',
    ]
    def __init__(self, document):
        QSyntaxHighlighter.__init__(self, document)
        tri = (quote)
        trid = (dquote)
        # Multi-line strings (expression, flag, style)
        # FIXME: The triple-quotes in these two lines will mess up the
        # syntax highlighting from this point onward
        self.tri_single = (QRegExp(tri), 1, STYLES['string2'])
        self.tri_double = (QRegExp(trid), 2, STYLES['string2'])

        rules = []

        # Keyword, operator, and brace rules
        rules += [(r'\b%s\b' % w, 0, STYLES['keyword'])
            for w in Highlighter.keywords]
        rules += [(r'%s' % o, 0, STYLES['operator'])
            for o in Highlighter.operators]
        rules += [(r'%s' % b, 0, STYLES['brace'])
            for b in Highlighter.braces]

        # All other rules
        rules += [
            # Numeric literals
            (r'\b[+-]?[0-9]+[lL]?\b', 0, STYLES['numbers']),
            (r'\b[+-]?0[xX][0-9A-Fa-f]+[lL]?\b', 0, STYLES['numbers']),
            (r'\b[+-]?[0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?\b', 0, STYLES['numbers']),

            # 'self'
            (r'\bself\b', 0, STYLES['self']),

            # Double-quoted string, possibly containing escape sequences ### "\"([^\"]*)\"" ### "\"(\\w)*\""
            (r'"[^"\\]*(\\.[^"\\]*)*"', 0, STYLES['string']),
            # Single-quoted string, possibly containing escape sequences
            (r"'[^'\\]*(\\.[^'\\]*)*'", 0, STYLES['string']),

            # 'def' followed by an word
            (r'\bdef\b\s*(\w+)', 1, STYLES['defclass']), ### (r'\bdef\b\s*(\w+)', 1, STYLES['defclass']),

            # 'self.' followed by an word
            #(r'\bself\b)', 1, STYLES['selfnext']), ### (r'\bself.\b\s*(\w+)', 1, STYLES['selfnext']), 

            # 'Q' followed by an word
            (r'\b[Q.]\b\s*(\w+)', 1, STYLES['Qnext']),

            # 'class' followed by an identifier
            (r'\bclass\b\s*(\w+)', 1, STYLES['classes']),

            # From '#' until a newline
            (r'#[^\n]*', 0, STYLES['comment']),

            # 'Q'  word
            #(r'\\bQ[A-Za-z]+\\b', 1, STYLES['Qtclass']), #(QRegExp("\\bQ[A-Za-z]+\\b")
        ]

        # Build a QRegExp for each pattern
        self.rules = [(QRegExp(pat), index, fmt)
            for (pat, index, fmt) in rules]

    def highlightBlock(self, text):
        # Apply syntax highlighting to the given block of text.

        # Do other syntax formatting
        for expression, nth, format in self.rules:
            match = expression.match(text)

            while match.hasMatch():
                # Get the index of the nth match
                index = match.capturedStart(nth)
                length = len(match.captured(nth))
                self.setFormat(index, length, format)
                match = expression.match(text, index + length)

        self.setCurrentBlockState(0)


    def highlightBlock_old(self, text):

#        Apply syntax highlighting to the given block of text.

        # Do other syntax formatting
        for expression, nth, format in self.rules:
            index = expression.indexIn(text, 0)

            while index >= 0:
                # We actually want the index of the nth match
                index = expression.pos(nth)
                length = len(expression.cap(nth))
                self.setFormat(index, length, format)
                index = expression.indexIn(text, index + length)

        self.setCurrentBlockState(0)

        # Do multi-line strings
        in_multiline = self.match_multiline(text, *self.tri_single)
        if not in_multiline:
            in_multiline = self.match_multiline(text, *self.tri_double)


    def match_multiline(self, text, delimiter, in_state, style):
        '''Do highlighting of multi-line strings. ``delimiter`` should be a
        ``QRegExp`` for triple-single-quotes or triple-double-quotes, and
        ``in_state`` should be a unique integer to represent the corresponding
        state changes when inside those strings. Returns True if we're still
        inside a multi-line string when this function is finished.
        '''
        # If inside triple-single quotes, start at 0
        if self.previousBlockState() == in_state:
            start = 0
            add = 0
        # Otherwise, look for the delimiter on this line
        else:
            start = delimiter.indexIn(text)
            # Move past this match
            add = delimiter.matchedLength()

        # As long as there's a delimiter match on this line...
        while start >= 0:
            # Look for the ending delimiter
            end = delimiter.indexIn(text, start + add)
            # Ending delimiter on this line?
            if end >= add:
                length = end - start + add + delimiter.matchedLength()
                self.setCurrentBlockState(0)
            # No; multi-line string
            else:
                self.setCurrentBlockState(in_state)
                length = len(text) - start + add
            # Apply formatting
            self.setFormat(start, length, style)
            # Look for the next match
            start = delimiter.indexIn(text, start + length)

        # Return True if still inside a multi-line string, False otherwise
        if self.currentBlockState() == in_state:
            return True
        else:
            return False
