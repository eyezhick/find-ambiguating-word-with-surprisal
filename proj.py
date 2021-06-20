import math
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import numpy as np

class Symbol():
    def __init__(self, sym, s1=-1, s2=-1):
        self.sym = sym
        self.s1 = s1
        self.s2 = s2

    def __eq__(self, other):
        return self.sym == other.sym

    def __str__(self):
        if self.s1 == -1:
            return self.sym
        return self.sym + "_" + str(self.s1) + "-" + str(self.s2)

    def __repr__(self):
        return self.__str__()

    def is_terminal(self):
        raise Exception("is_terminal not implemented")


class NT(Symbol):
    def __str__(self):
        if self.s1 == -1:
            return "{}".format(self.sym)
        return "{}{}".format(self.sym, "_" + str(self.s1) + "-" + str(self.s2))

    def is_terminal(self):
        return False


class T(Symbol):
    def __str__(self):
        return "{}".format(self.sym)

    def is_terminal(self):
        return True


class Rule():
    def __init__(self, lhs, v, rhs):
        self.lhs = lhs
        self.rhs = rhs
        self.v = v

    def __str__(self):
        return "{}, {} => {}".format(self.v, self.lhs, self.rhs)

    def __repr__(self):
        return self.__str__()

    def target(self):
        return self.lhs

    def output(self):
        return self.rhs

    def value(self):
        return self.v

    def is_terminal(self):
        return all([o.is_terminal() for o in self.output()])

    def __eq__(self, other):
        if str(self.target()) == str(other.target()) and str(self.rhs[0]) == str(other.output()[0]) \
                and str(self.rhs[1]) == str(other.output()[1]):
            return True
        return False


# assume no ambiguity
def find_terminal_rule(pcfg, word):
    rules = []
    for rule in pcfg:
        if rule.output()[0] == word:
            rules.append(rule)
    return rules


def intersect(pcfg, words):
    '''
    Intersect the given pcfg with the the given sequence of words.
    :param pcfg: original pcfg
    :param words: list of terminal symbols
    :return: weighted cfg
    '''
    initial_symbol = pcfg[0].target()
    new_pcfg = []
    NonTerminals = []
    length = len(words)
    # tag each words with T
    tagged_words = [T(word) for word in words]

    # set of terminal rules
    for i, word in enumerate(tagged_words):
        rules = find_terminal_rule(pcfg, word)
        for rule in rules:
            NonTerminals.append(NT(rule.target().sym, i, i+1))
            new_pcfg.append(Rule(NT(rule.target().sym, i, i+1),
                                 rule.value(),
                                 rule.output()))

    for rule in pcfg:
        if rule.is_terminal():
            NonTerminals.append(NT(rule.target().sym, length, length))
            new_pcfg.append(Rule(NT(rule.target().sym, length, length),
                                 rule.value(),
                                 rule.output()))

    clear = True
    while clear:
        clear = False
        for rule in pcfg:
            # iterate nonterminal
            if not rule.is_terminal():
                for nt1 in NonTerminals:
                    for nt2 in NonTerminals:

                        if rule.output()[0] == nt1 and rule.output()[1] == nt2 and nt1.s2 == nt2.s1:
                            new_rule = Rule(NT(rule.target().sym, nt1.s1, nt2.s2),
                                            rule.value(),
                                            [nt1, nt2])
                            #print(new_rule)
                            if new_rule in new_pcfg:
                                continue
                            new_pcfg.append(new_rule)
                            NonTerminals.append(NT(rule.target().sym, nt1.s1, nt2.s2))
                            clear = True

    #print(pretty_print(new_pcfg))
    return remove_unreachable_rules(new_pcfg, NT(initial_symbol.sym, 0, length))


def remove_unreachable_rules(pcfg, initial):
    '''
    Remove unreachable rules from the given pcfg.
    Also orders the rules correctly (initial symbol the the lhs of 1st rule).

    :param pcfg: pcfg
    :param initial: initial nonterminal symbol
    :return: pcfg with the unreachable rules removed and correctly ordered
    '''
    nont = [initial]
    new_pcfg = []
    for nt in nont:
        for rule in pcfg:
            if not rule.is_terminal():
                if str(nt) == str(rule.target()):
                    if rule not in new_pcfg:
                        new_pcfg.append(rule)
                        if not [nt for nt in nont if str(nt) == str(rule.output()[0])]:
                            nont.append(rule.output()[0])
                        if not [nt for nt in nont if str(nt) == str(rule.output()[1])]:
                            nont.append(rule.output()[1])
            else:
                if str(nt) == str(rule.target()):
                    if rule not in new_pcfg:
                        new_pcfg.append(rule)

    return new_pcfg


def inside(pcfg, symbol):
    '''
    Calculate inside value. As assumed in class, inside*(n) = 1.
    :param pcfg: the weighted cfg
    :param symbol: initial nonterminal symbol
    :return: total weight of the given weighted cfg
    '''
    if symbol.is_terminal():
        return 1
    if symbol.s1 == symbol.s2:
        return 1

    sum = 0
    x = 0
    for rule in pcfg:
        if str(rule.target()) == str(symbol):
            if len(rule.output()) == 1:
                # is terminal rule
                sum = sum + rule.value()
            else:
                if str(symbol) == str(rule.output()[0]):
                    # need some algebra here
                    x = x + rule.value() * inside(pcfg, rule.output()[1])
                else:
                    sum = sum + (rule.value() * inside(pcfg, rule.output()[0]) * inside(pcfg, rule.output()[1]))

    return sum/(1-x)


def surprisal(pcfg, words):
    '''
    Calculate surprisal. intersected grammars and their total weights are calcuated.

    :param pcfg: original pcfg
    :param words: list of terminal symbols
    :return: surprisal value
    '''
    prev = words[:-1]
    current = words

    g1_weighted = intersect(pcfg, prev)
    g1_total_weight = inside(g1_weighted, g1_weighted[0].target())

    g2_weighted = intersect(pcfg, current)
    g2_total_weight = inside(g2_weighted, g2_weighted[0].target())

    return math.log(g1_total_weight/g2_total_weight, 2)


def pretty_print(pcfg):
    '''
    print weighted cfg in the same way as it is presented in the homework

    :param pcfg: weighted pcfg
    :return: None
    '''
    for rule in pcfg:
        if rule.is_terminal():
            print(f"{round(rule.value(), 6): <15}", rule.target(), " --> ", "\"" + str(rule.output()[0]) + "\"")
        else:
            print(f"{round(rule.value(), 6): <15}", rule.target(), " --> ", rule.output()[0], rule.output()[1])
    print()


def pretty_print2(words, value):
    '''
    print the surprisal in a clean format
    :param words: list of terminal symbols
    :param value: surprisal value
    :return: None
    '''
    string = ' '.join(words)
    print(f"{string: <30}", "\nsurprisal: ", round(value, 6))


def pretty_print3(ambi, disambi):
    '''
    print the ambiguating words in a clean format
    :param ambi: list of tuple containing the words and their indices
    :param disambi: disambiguating word
    :return: None
    '''
    print("Disambiguating word: \n\t" + disambi)
    for word in ambi:
        print("List of ambiguating words")
        print("\t\"" + str(word[0]) + "\"" + " at index " + str(word[1]))


def generate_barplot(pcfg, string, title=""):
    fig = plt.figure()
    ax = fig.add_axes([0, 0, 1, 1])
    y = []
    for ind, word in enumerate(string, 1):
        y.append(surprisal(pcfg, string[:ind]))

    print(string)
    print(y)
    ax.bar(string, y)
    plt.title(title)
    plt.show()

    #plt.savefig("./mygraph.png")


def get_terminal_rules(pcfg):
    terminal_syms = []
    for rule in pcfg:
        if rule.is_terminal():
            terminal_syms.append(rule)
    return terminal_syms


def get_ambiguous_terminals(terminals):
    ambi = []
    for i1, t1 in enumerate(terminals):
        for i2 in range(i1, len(terminals)):
            t2 = terminals[i2]
            if (t1.target().s1 == t2.target().s1) and (t1.target().s2 == t2.target().s2) and \
                    (t1.target().sym == t2.target().sym) and (t1.target().s1 != t1.target().s2):
                ambi.append(t1)
    return ambi


def ambiguity_finder(pcfg, string):
    '''
    Pseudocode
        1.	Calculate the surprisal values for all prefixes.
        2.	Locate the word with highest surprisal prediction. This is the disambiguating word.
        3.	Create two weighted CFGs; one intersected with the string up to the disambiguating word,
            and an-other one intersected with the string up to just be-fore the disambiguating word.
        4.	Filter out unreachable rules in each CFGs.
        5.	Retrieve the terminal rules that output the same word at the same index with different target NT symbol.
            Repeat for both CFGs.
        6.	Isolate the rules unique to on CFG. The terminal symbols targeted by these rules are the ambiguat-ing words.


    :return: a list of tuple contain the terminal symbol and its index.
             disambiguating word (string).
    '''

    # find the disambiguating word
    surprisal_predictions = []
    for ind, word in enumerate(string, 1):
        surprisal_predictions.append(surprisal(pcfg, string[:ind]))
    disam_index = surprisal_predictions.index(max(surprisal_predictions))
    wcfg1 = intersect(pcfg, string[:disam_index])
    wcfg2 = intersect(pcfg, string[:disam_index+1])
    terminals1 = get_terminal_rules(wcfg1)
    terminals2 = get_terminal_rules(wcfg2)
    ambi_terminal1 = get_ambiguous_terminals(terminals1)
    ambi_terminal2 = get_ambiguous_terminals(terminals2)

    res = []
    for t1 in ambi_terminal1:
        check = True
        for t2 in ambi_terminal2:
            if str(t1) == str(t2):
                check = False
        if check:
            res.append((t1.output()[0], t1.target().s2))

    return res, string[disam_index]


#######################################################################################################################
####################################### NP/Z ambiguity: find ambiguating word #########################################
#######################################################################################################################
def main():

    NPZ = ["the", "banker", "told", "about", "the", "buy-back", "resigned"]

    print("NP/Z sentence: \n\t" + ' '.join(NPZ))
    ambi, disam = ambiguity_finder(pcfg_npz, NPZ)
    pretty_print3(ambi, disam)


######################################### Example pCFGs used in the papaer ############################################

#Section 2.1
the_banker1 = ["the", "banker", "told", "about", "the", "buy-back", "resigned"]
the_banker2 = ["the", "banker", "told", "about", "the", "boss", "who", "was",
               "told", "about", "the", "buy-back", "resigned"]

pcfg_re = [
        Rule(NT("S"), 0.574928, [NT("NP"), NT("VP")]),
        Rule(NT("S"), 0.425072*0.11043, [NT("VBD"), NT("PP")]),
        Rule(NT("S"), 0.425072*0.141104, [NT("VBD"), NT("NPPP")]),
        Rule(NT("S"), 0.425072*0.214724, [NT("AUX"), NT("VP")]),
        Rule(NT("S"), 0.425072*0.484663, [NT("VBN"), NT("PP")]),
        Rule(NT("S"), 0.425072*0.0490798*0.74309393, [T("told")]),
        Rule(NT("S"), 0.425072*0.0490798*0.25690607, [T("resigned")]),
        Rule(NT("SBAR"), 1.0, [NT("WHNP"), NT("VP")]),
        Rule(NT("NP"), 0.8041237, [NT("DT"), NT("NN")]),
        Rule(NT("NP"), 0.0824742, [NT("NP"), NT("SBAR")]),
        Rule(NT("NP"), 0.11340206, [NT("NP"), NT("VP")]),
        Rule(NT("VP"), 0.11043, [NT("VBD"), NT("PP")]),
        Rule(NT("VP"), 0.141104, [NT("VBD"), NT("NPPP")]),
        Rule(NT("NPPP"), 1.0, [NT("NP"), NT("PP")]),
        Rule(NT("VP"), 0.214724, [NT("AUX"), NT("VP")]),
        Rule(NT("VP"), 0.484663, [NT("VBN"), NT("PP")]),
        Rule(NT("VP"), 0.0490798*0.74309393, [T("told")]),
        Rule(NT("VP"), 0.0490798*0.25690607, [T("resigned")]),
        Rule(NT("PP"), 1.0, [NT("IN"), NT("NP")]),
        Rule(NT("WHNP"), 1.0, [T("who")]),
        Rule(NT("DT"), 1.0, [T("the")]),
        Rule(NT("NN"), 0.33, [T("boss")]),
        Rule(NT("NN"), 0.33, [T("banker")]),
        Rule(NT("NN"), 0.33, [T("buy-back")]),
        Rule(NT("IN"), 0.5, [T("about")]),
        Rule(NT("IN"), 0.5, [T("by")]),
        Rule(NT("AUX"), 1.0, [T("was")]),
        Rule(NT("VBD"), 0.74309393, [T("told")]),
        Rule(NT("VBD"), 0.25690607, [T("resigned")]),
        Rule(NT("VBN"), 1.0, [T("told")])

]

# subject/object relative cluase example. Section 2.1
the_man_who_saw_you = ["the", "man", "who", "saw", "you", "saw", "me"]
the_man_who_you_saw = ["the", "man", "who", "you", "saw", "saw", "me"]

pcfg3 = [
        Rule(NT("S"), 1.0, [NT("NP"), NT("VP")]),
        Rule(NT("S+R"), 0.5, [NT("NP+R"), NT("VP")]),
        Rule(NT("S+R"), 0.5, [NT("NP+R"), NT("S/NP")]),
        Rule(NT("NP"), 0.33, [NT("DT"), NT("NBAR")]),
        Rule(NT("NP"), 0.33, [T("you")]),
        Rule(NT("NP"), 0.33, [T("me")]),
        Rule(NT("NBAR"), 0.5, [NT("NBAR"), NT("S+R")]),
        Rule(NT("NBAR"), 0.5, [T("man")]),
        Rule(NT("S/NP"), 1.0, [NT("NP"), NT("V")]),
        Rule(NT("VP"), 1.0, [NT("V"), NT("NP")]),
        Rule(NT("V"), 0.33, [T("saw")]),
        Rule(NT("V"), 0.33, [T("told")]),
        Rule(NT("V"), 0.33, [T("passed")]),
        Rule(NT("NP+R"), 1.0, [T("who")]),
        Rule(NT("DT"), 1.0, [T("the")])

]

# NP/Z example. Section 3

pcfg_npz = [
        Rule(NT("S"), 0.574928+0.01550275902+0.00535968971, [NT("NP"), NT("VP")]),
        Rule(NT("S"), 0.425072*0.11043, [NT("VBD"), NT("PP")]),
        Rule(NT("S"), 0.425072*0.141104, [NT("VBD"), NT("NPPP")]),
        Rule(NT("S"), 0.425072*0.214724, [NT("AUX"), NT("VP")]),
        Rule(NT("S"), 0.425072*0.484663, [NT("VBN"), NT("PP")]),
        Rule(NT("SBAR"), 1.0, [NT("WHNP"), NT("VP")]),
        Rule(NT("NP"), 0.8041237, [NT("DT"), NT("NN")]),
        Rule(NT("NP"), 0.0824742, [NT("NP"), NT("SBAR")]),
        Rule(NT("NP"), 0.11340206, [NT("NP"), NT("VP/N")]),
        Rule(NT("VP"), 0.22, [NT("VBD"), NT("PP")]),
        Rule(NT("VP"), 0.25, [NT("VBD"), NT("NPPP")]),
        Rule(NT("NPPP"), 1.0, [NT("NP"), NT("PP")]),
        Rule(NT("VP"), 0.43, [NT("AUX"), NT("VP")]),
        Rule(NT("VP/N"), 1, [NT("VBN"), NT("PP")]),
        Rule(NT("VP"), 0.1*0.74309393, [T("told")]),
        Rule(NT("VP"), 0.1*0.25690607, [T("resigned")]),
        Rule(NT("PP"), 1.0, [NT("IN"), NT("NP")]),
        Rule(NT("WHNP"), 1.0, [T("who")]),
        Rule(NT("DT"), 1.0, [T("the")]),
        Rule(NT("NN"), 0.33, [T("boss")]),
        Rule(NT("NN"), 0.33, [T("banker")]),
        Rule(NT("NN"), 0.33, [T("buy-back")]),
        Rule(NT("IN"), 0.5, [T("about")]),
        Rule(NT("IN"), 0.5, [T("by")]),
        Rule(NT("AUX"), 1.0, [T("was")]),
        Rule(NT("VBD"), 0.74309393, [T("told")]),
        Rule(NT("VBD"), 0.25690607, [T("resigned")]),
        Rule(NT("VBN"), 1.0, [T("told")])

]

if __name__ == "__main__":
    main()
