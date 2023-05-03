// Percy++ Copyright 2007,2012,2013,2014
// Ian Goldberg <iang@cs.uwaterloo.ca>,
// Ann Yang <y242yang@uwaterloo.ca>,
// Casey Devet <cjdevet@uwaterloo.ca>,
// Paul Hendry <pshendry@uwaterloo.ca>,
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

#include <vector>
#include "subset.h"
#include "subset_iter.h"
#include "recover.h"

NTL_CLIENT

bool EasyRecover(dbsize_t bytes_per_word, nservers_t t,
	nservers_t h, vector<DecoderResult<GF2E> > &results,
	dbsize_t word_number, const vec_GF2E &values_vec, 
	const vec_GF2E &indices_vec, const vec_GF2E &interp_indices)
{
    //vector<DecoderResult<GF2E> > newresults;
    vector<DecoderResult<GF2E> >::const_iterator Riter;
    vector<vector<GF2E> > newF;
    vector<vector<nservers_t> > newG;

#ifdef VERBOSE_RECOVER
    std::cerr << "EasyRecover:\n";
#endif
    for (Riter = results.begin(); Riter != results.end(); ++Riter) {
#ifdef VERBOSE_RECOVER
    std::cerr << "  iteration:\n";
#endif
        // Pick a random subset I of G, of size t+1
        vector<nservers_t> I;
        unsigned short int deg = t;
	random_subset(Riter->G, I, deg+1);
        // Use Lagrange interpolation to find the unique polynomial phi
        // of degree t which matches the points indexed by I
        vec_GF2E I_indices, I_values;
        vector<nservers_t>::const_iterator Iiter;
        I_indices.SetLength(deg+1);
        I_values.SetLength(deg+1);
        nservers_t i = 0;
        for (Iiter = I.begin(); Iiter != I.end(); ++i, ++Iiter) {
            I_indices[i] = indices_vec[*Iiter];
            I_values[i] = values_vec[*Iiter];
        }
        GF2EX phi;
        phi = interpolate(I_indices, I_values);

        // Count the number of points in G that agree, and that
        // disagree, with phi
        nservers_t numagree = 0;
        nservers_t numdisagree = 0;
        newG.push_back(vector<nservers_t>());
        vector<nservers_t>::const_iterator Giter;
        for (Giter = Riter->G.begin(); Giter != Riter->G.end(); ++Giter) {
            if (eval(phi, indices_vec[*Giter]) == values_vec[*Giter]) {
                ++numagree;
                newG.back().push_back(*Giter);
            } else {
                ++numdisagree;
            }
        }

        // If at least h agreed, and less than h-t disagreed, then phi
        // can be the *only* polynomial of degree t matching at least
        // h points.
        if (numagree >= h && numdisagree < h-t) {
	    vector<GF2E> vec_wz;
	    GF2E wz;
	    if (interp_indices.length() == 0) {
		wz = eval(phi, GF2E::zero());
		vec_wz.push_back(wz);
	    } else {
		for (int i=0; i<interp_indices.length(); i++) {
		    wz = eval(phi, interp_indices[i]);
		    vec_wz.push_back(wz);
		}
	    }
            newF.push_back(vec_wz);
#ifdef VERBOSE_RECOVER
	    std::cerr << "        " << phi << "\n";
#endif
        } else {
            // This either isn't the right polynomial, or there may be
            // more than one.  Abort.
            return false;
        }
    }

    for (unsigned int i = 0; i < newF.size(); ++i) {
        results[i].G = newG[i];
	if (interp_indices.length() == 0) {
	    results[i].recovered.resize(1);
            results[i].recovered[0][word_number] = newF[i][0];
	} else {
	    results[i].recovered.resize(interp_indices.length());
	    for (unsigned int j = 0; j < interp_indices.length(); j++) {
		results[i].recovered[j][word_number] = newF[i][j];
	    }
	}
    }

    return true;
}


bool EasyRecover(dbsize_t bytes_per_word, nservers_t t,
	nservers_t h, vector<DecoderResult<ZZ_p> > &results, 
	dbsize_t word_number, const vec_ZZ_p &values,
	const vec_ZZ_p &indices, const vec_ZZ_p &interp_indices)
{
    vector<vector<ZZ_p> > newF;
    vector<vector<nservers_t> > newG;
    vector<DecoderResult<ZZ_p> >::const_iterator Riter;
    for (Riter = results.begin(); Riter != results.end();
            ++Riter) {
        // Pick a random subset I of G, of size t+1
        vector<nservers_t> I;

        nservers_t numagree, numdisagree;
        newG.push_back(vector<nservers_t>());

        ZZ_pX phi;
	random_subset(Riter->G, I, t+1);
	RSDecoder_ZZ_p::test_interpolate(t, values, indices, I, Riter->G,
	    numagree, numdisagree, newG.back(), phi);
        // If at least h agreed, and less than h-t disagreed, then phi
        // can be the *only* polynomial of degree t matching at least
        // h points.
        if (numagree >= h && numdisagree < h-t) {
            // Find the secret determined by phi
            vector<ZZ_p> vec_wz;
	    ZZ_p wz;
	    if (interp_indices.length() != 0){
	        for (int i=0; i<interp_indices.length(); i++) {
                    eval(wz, phi, interp_indices[i]);
	            vec_wz.push_back(wz);
	        }
	    } else {
		eval(wz, phi, ZZ_p::zero());
		vec_wz.push_back(wz);
	    }
            newF.push_back(vec_wz);
        } else {
            // This either isn't the right polynomial, or there may be
            // more than one.  Abort, and we'll use HardRecover.
            return false;
        }
    }
    for (unsigned int i = 0; i < newF.size(); ++i) {
        results[i].G = newG[i];
        if (interp_indices.length() != 0){
	    results[i].recovered.resize(interp_indices.length());
	    for (unsigned int j = 0; j < interp_indices.length(); j++) {
                results[i].recovered[j][word_number] = newF[i][j];
	    }
	} else {
	    results[i].recovered.resize(1);
	    results[i].recovered[0][word_number] = newF[i][0];
	}
    }

    return true;
}

