# Copyright (C) 2012, Joao Rodrigues (anaryin@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
Assorted utility functions to manipulate structural data.

This module serves as a container for all structural manipulation
utilities that used to belong in the structure objects themselves.

Functions should take as argument a collection of objects, or those
that can be unfolded in children (Structure, Chain, ...)
The output should be an equivalent data structure.

The module should also be callable directly to allow straightforward
manipulation in a shell enviroment.

Available functions

+ renumber_residues
"""

from Bio.PDB.Selection import unfold_entities

# Private renumbering method
#
def _renumber_residue(container, old_id, new_id):
    """
    To allow custom renumbering of residues.
    """
    # Identify residue
    residue = container[old_id]

    # Change old references in both dict and list
    del container.child_dict[residue.id]
    residue.id = new_id
    container.child_dict[residue.id] = residue

def renumber_residues(container, nstart=1, renum_method='conservative', chain_gap=0):
  """
  Changes the numbering of a collection of residues.
  
  Input: container (list of Residues, Chain, Structure, ...)
  Output: equivalent object (list of Residues, Chain, Structure, ...)
  
  Options:
    nstart (int)        - number to start renumbering from [default: 1]

    renum_method (str)  - Defines the (re)numbering mode to use:                         
                          + consecutive: simple consecutive renumbering from 'nstart'. Ignores
                            original numbering gaps.
                          + conservative: subtracts a seed from the residue number to make the
                            numbering start from 'nstart'. Preserves gaps. [DEFAULT]

    chain_gap (int)     - Numbering gap to introduce when a new chain is found
  """

  def _consecutive_renumbering(reslist, seed, chain_gap):
    """
    Renumbers a residue list consecutively
    starting from a predetermined seed.
    Ignores any gaps present in the original sequence.
    """   

    cur_chain = reslist[0].parent.id
    rindex = seed
    for res in reslist:
      if cur_chain != res.parent.id and chain_gap:
        cur_chain = res.parent.id
        rindex += chain_gap
      _renumber_residue(res.parent, res.id, (res.id[0], rindex, res.id[2]))      
      rindex += 1

  def _conservative_renumbering(reslist, iseed, chain_gap):
    """
    Renumbers a residue list consecutively
    starting from a predetermined seed.
    Keeps gaps present in the original sequence.
    """   

    i_fres = reslist[0].id[1]
    cur_chain = reslist[0].parent.id
    seed = iseed-i_fres
    for res in reslist:
      if cur_chain != res.parent.id:
        cur_chain = res.parent.id
        i_fres = res.id[1] 
        seed = iseed-i_fres
        if chain_gap:
          seed += chain_gap

      _renumber_residue(res.parent, res.id, (res.id[0], res.id[1]+seed, res.id[2]))

  # Check input
  allowed_types = ['S','M','I','C']
  try:
    assert container.level in allowed_types
    residues = unfold_entities(container, 'R')
  except AttributeError: # Otherwise, list of residues?
    try: 
      for element in container:
        assert element.level == "R"
      residues = container
    except (TypeError, AttributeError, AssertionError):
      msg = "Input has to be either a Biopython Structure, Model, Interface, "
      msg+= "or Chain object, or otherwise an iterable of Residue objects."
      raise Exception(msg)
    except Exception, e:
      msg = "An unexpected error occured: %s" %e
      raise Exception(msg)

  # Find & apply renumbering method
  try:
    retcode = locals()['_%s_renumbering' %renum_method](residues, nstart, chain_gap)
  except KeyError:
    msg = "Numbering mode unrecognized: %s\n" %renum_method
    msg+= "Use 'consecutive' or 'conservative'.\n"
    msg+= "Check the documentation for further information."
    raise Exception(msg)
  return retcode

if __name__ == '__main__':
  pass
     