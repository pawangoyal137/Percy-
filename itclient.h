// Percy++ Copyright 2007,2012,2013,2014
// Ian Goldberg <iang@cs.uwaterloo.ca>,
// Casey Devet <cjdevet@uwaterloo.ca>,
// Ann Yang <y242yang@uwaterloo.ca>
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

#ifndef __ITCLIENT_H__
#define __ITCLIENT_H__

#include <vec_vec_ZZ_p.h>
#include <vec_GF2E.h>
#include "percyclient.h"
#include "gf2e.h"
#include "percyparams.h"
#include "itparams.h"
#include "rsdecoder.h"

#ifdef SPIR_SUPPORT
#include "spirclient.h"
#endif

NTL_CLIENT

/// A PIR client for the IT-PIR protocol by Goldberg (2007) over the integers
/// modulo p.
/// This protocol was introduced in
/// <a href="http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=4223220">Improving the Robustness of Private Information Retrieval</a>
/// and revisited in
/// <ul>
/// <li><a href="http://link.springer.com/chapter/10.1007%2F978-3-642-27576-0_13">Revisiting the Computational Practicality of Private Information Retrieval</a></li>
/// <li><a href="http://dl.acm.org/citation.cfm?id=2046784">Practical PIR for Electronic Commerce</a></li>
/// <li><a href="https://www.usenix.org/conference/usenixsecurity12/technical-sessions/presentation/devet">Optimally Robust Private Information Retrieval</a></li>
/// <li><a href="http://internetsociety.org/doc/one-block-size-fits-all-pir-and-spir-variable-length-records-multi-block-queries">One (Block) Size Fits All: PIR and SPIR Over Arbitrary-Length Records via Multi-block PIR Queries</a></li>
/// <li><a href="http://dl.acm.org/citation.cfm?id=2517840.2517854">Outsourced Private Information Retrieval with Pricing and Access Control</a></li>
/// </ul>

class PercyClient_ZZ_p : public PercyClient {
public:
    /// Constructor.
    /// @param params	    Parameters for the client.
    /// @param num_servers  The number of servers used
    /// @param t	    The privacy level.  I.e. the maximum number of 
    ///			    servers that can collude and the queries will remain
    ///			    private.
    /// @param sids	    An array of IDs for the servers that will be used.
    /// @param stats	    Statistics collection object.  No statistics will be
    ///			    collected if NULL.
    PercyClient_ZZ_p (const PercyClientParams * params, nservers_t num_servers,
	    nservers_t t, sid_t * sids, PercyStats * stats = NULL);

    /// Destructor.
    virtual ~PercyClient_ZZ_p ();

private:
    // Virtual members as described in PercyClient class
    virtual void encode_request_impl (nqueries_t request_identifier);
    virtual dbsize_t send_request_impl (nqueries_t request_identifier, 
	    vector<ostream*> &osvec, bool send_num_queries = true);
    virtual dbsize_t receive_replies_impl (nqueries_t request_identifier,
	    vector<istream*> &isvec);
    virtual nqueries_t process_replies_impl (nservers_t h,
	    vector<vector<PercyResult> >& results);

    // Based on the sids of the servers, choose the index that will
    // correspond with each server.  This function is called by the
    // constructor.
    virtual void choose_indices(sid_t *sids);
    virtual void choose_interp_indices(nqueries_t qbs);

    // Private ZZ_p members
    const ZZ_pParams * params;
    vec_ZZ_p indices;
    vec_ZZ_p interp_indices;
    vector<vec_ZZ_p> vecs_interp_indices;
    map<nqueries_t, vec_vec_ZZ_p *> queries;
#ifdef SPIR_SUPPORT
    map<nqueries_t, SPIRClientQuery *> spir_query_info;
#endif
    map<nqueries_t, vec_ZZ_p> randmults;
    vector<vector<vec_ZZ_p> > answers;
    vector<vector<DecoderResult<ZZ_p> > > unfinished_results;
    vector<std::set<dbsize_t> > decoded;
};


/// A PIR client for the IT-PIR protocol by Chor et al.\ (1995).
/// This protocol was introduced in
/// <a href="http://dl.acm.org/citation.cfm?id=796270">Private Information Retrieval</a>.

class PercyClient_Chor : public PercyClient {
public:
    /// Constructor.
    /// @param params	    Parameters for the client.
    /// @param num_servers  The number of servers used
    /// @param stats	    Statistics collection object.  No statistics will be
    ///			    collected if NULL.
    PercyClient_Chor (const PercyClientParams * params, nservers_t num_servers,
	    PercyStats * stats = NULL);

    /// Destructor.
    virtual ~PercyClient_Chor ();

private:
    // Virtual members as described in PercyClient class
    virtual void encode_request_impl (nqueries_t request_identifier);
    virtual dbsize_t send_request_impl (nqueries_t request_identifier, 
	    vector<ostream*> &osvec, bool send_num_queries = true);
    virtual dbsize_t receive_replies_impl (nqueries_t request_identifier,
	    vector<istream*> &isvec);
    virtual nqueries_t process_replies_impl (nservers_t h,
	    vector<vector<PercyResult> >& results);

    // Private Chor members
    const ChorParams * params;
    map<nqueries_t, unsigned char *> queries;
    vector<unsigned char *> answers;
};


/// A PIR client for the IT-PIR protocol by Goldberg (2007) over GF(2^E).
/// Supported fields are GF(2^8) and GF(2^16).
/// This protocol was introduced in
/// <a href="http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=4223220">Improving the Robustness of Private Information Retrieval</a>
/// and revisited in
/// <ul>
/// <li><a href="http://link.springer.com/chapter/10.1007%2F978-3-642-27576-0_13">Revisiting the Computational Practicality of Private Information Retrieval</a></li>
/// <li><a href="http://dl.acm.org/citation.cfm?id=2046784">Practical PIR for Electronic Commerce</a></li>
/// <li><a href="https://www.usenix.org/conference/usenixsecurity12/technical-sessions/presentation/devet">Optimally Robust Private Information Retrieval</a></li>
/// <li><a href="http://internetsociety.org/doc/one-block-size-fits-all-pir-and-spir-variable-length-records-multi-block-queries">One (Block) Size Fits All: PIR and SPIR Over Arbitrary-Length Records via Multi-block PIR Queries</a></li>
/// <li><a href="http://dl.acm.org/citation.cfm?id=2517840.2517854">Outsourced Private Information Retrieval with Pricing and Access Control</a></li>
/// </ul>

template<typename GF2E_Element>
class PercyClient_GF2E : public PercyClient {
public:
    /// Constructor.
    /// @param params	    Parameters for the client.
    /// @param num_servers  The number of servers used
    /// @param t	    The privacy level.  I.e. the maximum number of 
    ///			    servers that can collude and the queries will remain
    ///			    private.
    /// @param sids	    An array of IDs for the servers that will be used.
    /// @param stats	    Statistics collection object.  No statistics will be
    ///			    collected if NULL.
    PercyClient_GF2E (const PercyClientParams * params, nservers_t num_servers,
	    nservers_t t, sid_t * sids, PercyStats * stats = NULL);

    /// Destructor.
    virtual ~PercyClient_GF2E ();

private:
    // Virtual members as described in PercyClient class
    virtual void encode_request_impl (nqueries_t request_identifier);
    virtual dbsize_t send_request_impl (nqueries_t request_identifier, 
	    vector<ostream*> &osvec, bool send_num_queries = true);
    virtual dbsize_t receive_replies_impl (nqueries_t request_identifier,
	    vector<istream*> &isvec);
    virtual nqueries_t process_replies_impl (nservers_t h,
	    vector<vector<PercyResult> >& results);

    virtual void choose_indices(sid_t *sids);
    virtual void choose_interp_indices(nqueries_t qbs);

    // A NTL-less method to attempt a fast recovery
    bool try_fast_recover (nservers_t h, vector<vector<PercyResult> > &results);
    // Some helpers for it
    void construct_lagrange_coeffs(GF2E_Element *coeffs, GF2E_Element alpha,
	    nservers_t firstpoint, nservers_t numpoints);
    inline GF2E_Element interpolate(const GF2E_Element *word_answers, 
	    const GF2E_Element *coeffs, nservers_t firstpoint,
	    nservers_t numpoints);

    // Private GF2E members
    const GF2EParams * params;
    GF2E_Element * indices;
    GF2E_Element * interp_indices;
    vec_GF2E vec_interp_indices;
    vector<vec_GF2E> vecs_interp_indices;
    vec_GF2E indices_ntl;
    map<nqueries_t, vector<GF2E_Element> > randmults;
    map<nqueries_t, GF2E_Element *> stored_shares;
    vector<GF2E_Element *> answers;
    vector<vector<vec_GF2E> > answers_ntl;
    map<nqueries_t, vector<dbsize_t> > undecoded_indices;
    vector<vector<DecoderResult<GF2E> > > unfinished_results;
    vector<std::set<dbsize_t> > decoded;
};

#include "itclient_impl.h"

#endif
