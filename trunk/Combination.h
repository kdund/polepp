#ifndef COMBINATION_H
#define COMBINATION_H

/**
   utility class to select combination of elements in a collection

   @author Tadashi Maeno
*/

#include <vector>

namespace AnalysisUtils {

  /** combination
   */
  template <class COLL> class Combination {

  public:

    /** constructor
	@param coll [in] collection
	@param nElement [in] number of element to be selected
     */
    Combination(COLL *coll, const unsigned int nElement);
    
    /** destructor
     */
    ~Combination() {}

    /** reset internal indices
     */
    void reset();

    /** get a combination. internal indices are incremented
	@param comb [out] combination of elements
	@retrun true:success false:end of selection
    */
    template <class OUT>
    bool get(OUT &comb)
    {
      // init returned combination
      comb.clear();
      
      if (m_first) {
	if (m_index[m_nElement-1] >= m_coll->size()) return false;
	m_first = false;
      } else {
	// increment indices
	if (!setNewIndex(m_nElement-1)) return false;
      }
      
      // assign to returned combination
      for (unsigned int i=0; i<m_nElement; ++i)
	comb.push_back((*m_coll)[m_index[i]]);
      
      return true;
    }

    /** get a combination which passes a selection criteria
	@param caller
	@param comb [out] combination of elements
	@param criteria [in] a function pointer of selection criteria 
	@retrun true:success false:end of selection
    */
    template <class CALLER, class OUT, class CRITERIA>
    bool goodOnes(CALLER *caller, OUT &comb, CRITERIA criteria)
    {
      bool ret = get(comb);
      if (!ret) return false;

      // check if this passes the criteria
      if (criteria(caller, comb)) return true;

      // if not, look for next combination
      return goodOnes(caller, comb, criteria);
    }

    /** get a combination and the remaining elements
     */
    template <class OUT>
    bool get(OUT &comb, OUT &remain)
    {
      // init returned vector
      remain.clear();

      bool ret = get(comb);
      if (!ret) return false;

      // assigin to remain
      unsigned int sIndex = 0;
      std::vector<unsigned int>::const_iterator it  = m_index.begin();
      std::vector<unsigned int>::const_iterator itE = m_index.end();
      for (; ; ++it)
	{
	  unsigned int eIndex;
	  if (it!=itE)
	    eIndex = *it;
	  else
	    eIndex = m_coll->size();

	  for (unsigned int i=sIndex; i<eIndex; ++i)
	    remain.push_back((*m_coll)[i]);
	  sIndex = eIndex+1;

	  if (it == itE) break;
	}

      return true;
    }

    /** get a combination and the remaining elements.
	the combination passes a selection criteria
     */
    template <class CALLER, class OUT, class CRITERIA>
    bool goodOnes(CALLER *caller, OUT &comb, OUT &remain, CRITERIA criteria)
    {
      bool ret = get(comb, remain);
      if (!ret) return false;

      // check if this passes the criteria
      if (criteria(caller, comb)) return true;

      // if not, look for next combination
      return goodOnes(caller, comb, remain, criteria);
    }

  private:

    //! collection
    COLL *m_coll;
    
    //! number of elements to be selected
    const unsigned int m_nElement;

    //! indices of elements
    std::vector<unsigned int> m_index;

    //! flag to check if first
    bool m_first;

    //! set new index recursively
    bool setNewIndex (int iElement);
  };


  ////////////////////////////////////////////////////
  //  implementation

  template <class COLL> inline Combination<COLL>::Combination(COLL *coll, const unsigned int nElement)
    : m_coll(coll), m_nElement(nElement), m_first(true)
  {
    // init indices
    //    internal indices are [0,1,2...] at beginning 
    for (unsigned int i=0; i<nElement; ++i)
      m_index.push_back(i);
  }
  
  template <class COLL> inline void Combination<COLL>::reset()
  {
    // reset indices
    m_index.clear();
    for (unsigned int i=0; i<m_nElement; ++i)
      m_index.push_back(i);

    // set first flag
    m_first = true;
  }

  // iElement runs over 0..(m_nElement-1)
  template <class COLL> inline bool Combination<COLL>::setNewIndex (int iElement)
  {
    if (iElement < 0) return false;
    
    if (iElement+1 == static_cast<int>(m_nElement))
      {
	if ((m_index[iElement]+1) < m_coll->size())
	  {
	    ++(m_index[iElement]);
	    return true;
	  }
      }
    else
      {
	if ((m_index[iElement]+1) < m_index[iElement+1])
	  {
	    ++(m_index[iElement]);
	    return true;
	  }
      }
    
    if (setNewIndex (iElement-1))
      {
	m_index[iElement] = m_index[iElement-1]+1;
	return true;
      }
    
    return false;
  }

}; // end of AnalysisUtils

#endif // END OF ANALYSISTOOLS_ANALYSISCOMBINATION_H
