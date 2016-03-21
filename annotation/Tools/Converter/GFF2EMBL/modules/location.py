#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

class Span(object):
    
    def __init__(self, start, end = None):
        self.complement = False
        self.start_modifier = None
        self.end_modifier = None
        self.complement = False
        self.start = None
        self.end = None
        self.join = ".."  # TODO: support other join operators!!
        self.accession = ""
        self.set_start(start)
        self.set_end(end)
    
    def __repr__(self):
        output = ""
        start = self.end if self.complement else self.start
        end = self.start if self.complement else self.end
        start_mod = self.end_modifier if self.complement else self.start_modifier
        end_mod = self.start_modifier if self.complement else self.end_modifier
        if start != None:
            output += start_mod if start_mod else ""
            output += "%s" % (start+1)
        if end:
            output += self.join
            output += end_mod if end_mod else ""
            output += "%s" % end
        if self.accession:
            output = "%s:%s" % (self.accession, output)
        if self.complement:
            output = "complement(%s)" % output
        return output
    
    def _set(self, data):
        value = None
        modifier = None
        if data == None:
            return value, modifier
        elif data == "single":
            value = ""
        elif type(data) == type(1):
            value = data
        elif data[0] in ["<",">"]:
            modifier = data[0]
            value = int(data[1:])
        else:
            value = int(data)
        
        return value, modifier
    
    def set_accession(self, value):
        self.accession = value
    
    def set_any(self):
        self.join = "."
    
    def set_between(self):
        self.join = "^"
        if self.end - self.start != 1:
            print "WARNING: Setting operator 'between' for non-adjoining start-end values"
    
    def set_complement(self):
        temp = self.end
        self.end = self.start
        self.start = temp
    
    def set_end(self, value):
        self.end, self.end_modifier = self._set(value)
        if self.start != None and self.end != None and self.end < self.start:
            self.complement = True
    
    def set_start(self, value):
        self.start, self.start_modifier = self._set(value)
        
        if self.start != None and self.end != None and self.end < self.start:
            self.complement = True

class Location(object):
    """
    The location descriptor can be one of the following: 
    (a) a single base number, e.g. "467"
    (b) a site between two indicated adjoining bases, e.g. "123^124"
    (c) a single base chosen from within a specified range of bases (not allowed for new
        entries), e.g. "102.110"
    (d) the base numbers delimiting a sequence span, e.g. "240..565"
    (e) a remote entry identifier followed by a local location descriptor
        (i.e., a-d), e.g. "J00194.1:100..202"
    
    
    multiple locations can be joined as "join(12..78,134..202)".
    locations on the complementary strand can be presented as "complement(34..126)"
    'fuzzy' ends can be represented as "<1..>888"
    """
    
    def __init__(self, *args, **kwargs):
        self.spans = []
        for arg in args:
            if hasattr(arg, "start"): # Deal with BioPython FeatureLocation basically
                span = Span(str(arg.start))
                if hasattr(arg, "end"):
                    span.set_end(str(arg.end))
                if hasattr(arg, "strand") and arg.strand < 0:
                    span.set_complement()
                if hasattr(arg, "ref") and hasattr(arg, "ref_db") and arg.ref and arg.ref_db:
                    span.accession = "%s:%s" % (arg.ref_db, arg.ref)
                self.spans += [span]
            elif not self.spans:
                self.spans += [Span(arg)]
            elif arg == "any":
                self.spans[-1].set_any()
            elif arg == "between":
                self.spans[-1].set_between()
            elif self.spans[-1].end != None:
                self.spans += [Span(arg)]
            else:
                self.spans[-1].set_end(arg)
    
    def __repr__(self):
        output = ""
        if not self.spans:
            return "."
        if len(self.spans) > 1:
            output = "join(%s)" % (",".join(map(str, self.spans)))
        else:
            output = "".join(map(str, self.spans))
        
        return output

if __name__ == '__main__':
    
    print " ----- Feature Location implementation for EMBL format ----- "
    print
    print "Location(\"<23\", \"700\")               =>", Location("<23", "700")
    print "Location(\"<23\", 700, \">332\", \"<931\") =>", Location("<23", 700, ">332", "<931")
    print "Location(\">2321\", 700, \"332\")        =>", Location(">2321", 700, "332")
    print "Location(\">2321\", \"single\", \"332\")   =>", Location(">2321", "single", "332")
    print "Location(\">23\", 166, \"any\")          =>", Location(">23", 166, "any")
    print "Location(\"165\", 166, \"between\")      =>", Location("165", 166, "between")
    
    try:
        from Bio.SeqFeature import FeatureLocation, BeforePosition, AfterPosition
        print "location = FeatureLocation(BeforePosition(5), AfterPosition(33))"
        location = FeatureLocation(BeforePosition(5), AfterPosition(33))
        print "Location( location, 100, \"<50\")      =>", Location( location, 100, "<50")
    except:
        pass
    