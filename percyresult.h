// Percy++ Copyright 2007,2012,2013,2014
// Ian Goldberg <iang@cs.uwaterloo.ca>,
// Casey Devet <cjdevet@uwaterloo.ca>,
// Ryan Henry <rhenry@cs.uwaterloo.ca>
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of version 2 of the GNU General Public License as
// published by the Free Software Foundation.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// There is a copy of the GNU General Public License in the COPYING file
// packaged with this plugin; if you cannot find it, write to the Free
// Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301 USA

#ifndef __PERCYRESULT_H__
#define __PERCYRESULT_H__

#include <vector>
#include <string>
#include <iostream>
#include "percytypes.h"

using namespace std;

struct PercyResult {
    PercyResult(vector<nservers_t> G, string sigma) : G(G), sigma(sigma) {}
    vector<nservers_t> G;
    string sigma;

    void dump(ostream &out, unsigned int depth = 0) const {
	string leader(2*depth, ' ');
	out << leader << "PercyResult {" << endl;
	out << leader << "  G = <";
	for (vector<nservers_t>::const_iterator i=G.begin(); i != G.end(); ++i) {
	    out << " " << *i;
	}
	out << " >" << endl;
	out << leader << "  sigma = \"" << sigma << "\"" << endl;
	out << leader << "}" << endl;
    }
};

inline ostream& operator<<(ostream &out, const PercyResult &r) { r.dump(out); return out; }
inline ostream& operator<<(ostream &out, const vector<PercyResult> &r) {
    out << "<" << endl;
    for (vector<PercyResult>::const_iterator i=r.begin();
	    i != r.end(); ++i) {
	i->dump(out,1);
    }
    out << ">" << endl;
    return out;
}


struct PercyBlockResults {
    dbsize_t block_number;
    vector<PercyResult> results;

    void dump(ostream &out, unsigned int depth = 0) const {
	string leader(2*depth, ' ');
	out << leader << "PercyBlockResults {" << endl;
	out << leader << "  block_number = " << block_number << endl;
	out << leader << "  results = <" << endl;
	for (vector<PercyResult>::const_iterator i=results.begin();
		i != results.end(); ++i) {
	    i->dump(out, depth+2);
	}
	out << leader << "  >" << endl;
	out << leader << "}" << endl;
    }
};

inline ostream& operator<<(ostream &out, const PercyBlockResults &r) { r.dump(out); return out; }
inline ostream& operator<<(ostream &out, const vector<PercyBlockResults> &r) {
    out << "<" << endl;
    for (vector<PercyBlockResults>::const_iterator i=r.begin();
	    i != r.end(); ++i) {
	i->dump(out,1);
    }
    out << ">" << endl;
    return out;
}

#endif
