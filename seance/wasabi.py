from commands import getstatusoutput as run
from sys import stderr, exit

import webbrowser
import json


def wasabi(filename, name, url, user) :
    wid,wfile = _upload_to_wasabi(filename, name, url, user)

    webbrowser.open_new_tab(_build_url(url, wid, wfile))

def _build_url(url, wid, wfile) :
    return "%s?share=%s&file=%s" % (url, wid, wfile)

def _upload_to_wasabi(filename, name, url, user) :
    curl_command = 'curl --silent \
                    -F "action=save" \
                    -F "writemode=new" \
                    -F "name=%s" \
                    -F "userid=%s" \
                    -F "file=@%s" %s' % (name, user, filename, url)
    
    status,output = run(curl_command)

    if status != 0 :
        print >> stderr, "ERROR wasabi upload failed (curl err %d)" % status
        exit(1)

    tmp = json.loads(output)
    return tmp['id'], tmp['name']

