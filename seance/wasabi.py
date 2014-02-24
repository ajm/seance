from commands import getstatusoutput as run
from sys import stderr, exit

import webbrowser
import json


def wasabi(filename, name, url) :
    wid,wfile = _upload_to_wasabi(filename, name, url)

    webbrowser.open_new_tab(_build_url(url, wid, wfile))

def _build_url(url, wid, wfile) :
    return "%s?share=%s&file=%s" % (url, wid, wfile)

def _upload_to_wasabi(filename, name, url) :
    curl_command = 'curl --silent \
                    -F "action=save" \
                    -F "writemode=new" \
                    -F "name=%s" \
                    -F "userid=mbQdOD" \
                    -F "file=@%s" %s' % (name, filename, url)
    
    status,output = run(curl_command)

    if status != 0 :
        raise Exception("curl failed")

    #print curl_command
    #print output

    tmp = json.loads(output)
    return tmp['id'], tmp['name']

