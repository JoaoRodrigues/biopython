# Copyright 2010 by Joao Rodrigues.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
This module provides code to work with the WWW version of WHATIF
servers provided at http://swift.cmbi.kun.nl/.

Functions:
add_hydrogens        Add Hydrogens to a protein using the HTOPO service.
"""

import urllib
import os, tempfile

# Utility functions
def smcra_to_pdb(smcra_object, temp_dir='/tmp/'):
    """
    Since all servers work with PDB files, we can't use the SMCRA structure directly.
    Using PDBIO & tempfile to make it work seamlessly.
    """
    
    temp_path = tempfile.mktemp( '.pdb', dir=temp_dir )
    
    io = PDBIO()
    io.set_structure(smcra_object)
    io.save(temp_path)
    
    f = open(temp_path, 'r')
    pdb_data = f.read()
    f.close()
    
    os.remove(temp_path)
    
    return pdb_data

def pdb_to_smcra(pdb_contents, temp_dir='/tmp/'):
    
    temp_path = tempfile.mktemp( '.pdb', dir=temp_dir )
    f = open(temp_path, 'r')
    f.write(pdb_contents)
    f.close()
        
    p = PDBParser()
    smcra_object = p.get_structure(temp_path)
 
    os.remove(temp_path)
    
    return smcra_object

# Dummy REST webservice to Test WHATIF servers

def whatif_test(): # Check if service is up
    pass

# HTOPO Service to Add Hydrogens

def add_hydrogens(data):
    
    if not isinstance(data, str): # Hope it's a SMCRA object..
        data = smcra_to_pdb(data)
        smcra_data = 1
    else:
        smcra_data = 0
    
    h_added = whatif(data, "htopo")
    
    if smcra_data:
        return pdb_to_smcra(h_added)
    else:
        return h_added

# Waiting for G.Vriend to add REST access
# Meanwhile using Sjoerd de Vries code. Copyright as follows:
## Copyright 2008, 2009 Sjoerd de Vries
## This file is part of the Spyder module: "bio", "atom"
## For licensing information, see LICENSE.txt
## Sjoerd de Vries, December 2006, ADAPTED November 2007 for the HADDOCK server

import httplib, mimetypes, urlparse

site = 'http://swift.cmbi.ru.nl/'
whatifcgi = 'wiw-cgi/FiledownCGI.py'

def whatif(pdbdata, request, outputname='hadded.pdb'):
  url = site + whatifcgi
  
  urlparts = urlparse.urlsplit(url)
  host = urlparts[1]
  selector = urlparts[2]
    
  try:
    loc = post_multipart(pdbdata, request)
    fields = []
    for l in loc.splitlines():
      l = l.strip()
      if l.startswith("<INPUT TYPE"):
        name = l[l.index("NAME=")+len("NAME="):].strip()[1:]
        name = name[:name.index('"')]
        val = l[l.index("VALUE=")+len("VALUE="):].strip()[1:]
        val = val[:val.index('"')]
        fields.append((name,val))
    for f in fields:
      if f[0] == 'ID':
        location = "/servers/tmp/" + f[1]
        break
    
  except Exception, e:
    raise Exception("Could not submit job to the WHATIF server\nReason: %s" %e)
  
  try:  
    return urllib.urlopen(site+ location + "/" + outputname).read()
  except:
    raise Exception("Cannot download WHATIF file")

# Misc Functions to build the request

def post_multipart(filevalue, request):
    url = site + whatifcgi
    fields = (('request', request), ('&PDB1', ''))
    filekey = '&FIL1'
    
    urlparts = urlparse.urlsplit(url)
    host = urlparts[1]
    selector = urlparts[2]
    
    fil = (filekey, "dummy", filevalue)
    
    content_type, body = encode_multipart_formdata(fields, [fil])
    counter = 0
    while 1:
      counter += 1
      ok = True
      try:
        h = httplib.HTTPConnection(host)
        headers = {
          'User-Agent': 'anonymous',
          'Content-Type': content_type    
        }
        h.request('POST', selector, body, headers)
        res = h.getresponse() 
      except:
        ok = False
        if counter == 5: raise 
      if ok: break
    return res.read()

def encode_multipart_formdata(fields, files):
    """
    Based on some code snippet found on the web in 2006 or so
    """
    BOUNDARY = '----------1234567890abcdefghij_$'
    CRLF = '\r\n'
    L = []
    for file in files:
      (filekey, filename, filevalue) = file
      L.append('--' + BOUNDARY)
      L.append('Content-Disposition: form-data; name="%s"; filename="%s"' % (filekey, filename))
      L.append('Content-Type: text/plain')
      L.append('')
      for line in filevalue.split('\n'):
        L.append(line)
    for (key, value) in fields:
        if value == None: value = ""
        L.append('--' + BOUNDARY)
        L.append('Content-Disposition: form-data; name="%s"' % str(key))
        L.append('')
        L.append(str(value))
    L.append('--' + BOUNDARY + '--')
    L.append('')
    body = CRLF.join(L)
    content_type = 'multipart/form-data; boundary=%s' % BOUNDARY
    return content_type, body


if __name__=="__main__":
    import sys
    from Bio.PDB import PDBIO, PDBParser

    if len(sys.argv) != 3:
        print "Expects 2 arguments:"
        print "\t - Input PDB File"
        print "\t - Output PDB File Name"

        sys.exit()
    
    f = open(sys.argv[1], 'r')
    pdb_data = f.read()
    f.close()
    
    print "Adding hydrogen atoms..."
    h_add = add_hydrogens(pdb_data)
    print "Success!"
    
    out = open(sys.argv[2], 'w')
    out.write(h_add)
    out.close()
    
    