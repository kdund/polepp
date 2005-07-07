#ifndef ANALYSISUTILS_ANALYSISPERMUTATION_H
#define ANALYSISUTILS_ANALYSISPERMUTATION_H

/**
   utility class to return permutation of a collection

   @author Tadashi Maeno
*/

#include <algorithm>
#include <vector>

#include "Combination.h"

namespace AnalysisUtils {

  /** permutation
   */
  template <class COLL> class Permutation
  {
  public:

    /** constructor
	@param coll collection
     */
    Permutation(COLL *coll, const unsigned int nElement)
      : m_coll(coll), m_first(true)
    {
      // init indices
      for (unsigned int i=0; i<coll->size(); ++i)
	m_index_for_comb.push_back(new unsigned int(i));

      m_comb = new Combination<std::vector<unsigned int *> > (&m_index_for_comb, nElement);
    }
    
    /** destructor
     */
    ~Permutation()
    {
      delete m_comb;

      for (unsigned int i=0; i<m_index_for_comb.size(); ++i)
	delete m_index_for_comb[i];
    }

    /** get a permutation. This method changes the sequence in increasing order
	and return true if succeeds. If the current sequence is the last permutation,
	return false
	@param perm a vector for permutation
     */
    template <class OUT>
    bool get(OUT &perm)
    {
      // init returned vector
      perm.clear();

      // get combination if first
      if (m_first)
	{
	  m_first=false;
	  if (! m_comb->get(m_index)) return false;
	  std::sort(m_index.begin(), m_index.end());
	}
      else
	{
	  // change sequence. if this is the last permutation. get next combination
	  if (! std::next_permutation (m_index.begin(), m_index.end()))
	    {
	      if (! m_comb->get(m_index)) return false;
	      std::sort(m_index.begin(), m_index.end());
	    }
	}
      
      // assign
      std::vector<unsigned int *>::const_iterator it  = m_index.begin();
      std::vector<unsigned int *>::const_iterator itE = m_index.end();
      for (; it!=itE; ++it)
	perm.push_back((*m_coll)[**it]);

      return true;
    }

    /** get a permutation which passes a selection criteria 
	@param perm a vector for permutation
	@param criteria a function pointer of selection criteria
     */
    template <class CALLER, class OUT, class CRITERIA>
    bool goodOnes(CALLER *caller, OUT &perm, CRITERIA criteria)
    {
      get(perm);
      
      // check if this passes the criteria
      if (criteria(caller,perm)) return true;

      // if not, look for next combination
      return goodOnes(caller, perm, criteria);
    }
      

  private:

    //! combination
    Combination<std::vector<unsigned int *> > *m_comb;

    //! collection
    COLL *m_coll;

    //! indices of elements
    std::vector<unsigned int *> m_index_for_comb;
    std::vector<unsigned int *> m_index;

    //! flag to check if first
    bool m_first;

  };
}; // end of AnalysisUtils

#endif
