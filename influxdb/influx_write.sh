#! /usr/bin/env bash

if [ "$#" -lt 1 ]; then
    echo "usage: $0 [data]"
    exit 1
fi

#eicdesk04
#TOKEN="FNDVnu7yKPbw6zEZ744Unkm3___bEhEnrEC95OsDckhtKZVdb6MjivjZxIU-hTlTEox6PPXMU2PQ0bGdpTRrag=="

#eicdesk03
TOKEN="0jdyWS61SWPxUVfzG-s7jxEoeyvFKc8nkfeGhK_9ztRrb5mn5J7wxUpW2P8kSHOJNg35pkbz_K4e9ldFnv0tBA=="

DATA="$@"

curl --request POST \
     "http://localhost:8086/api/v2/write?org=epic&bucket=epic&precision=ms" \
     --header "Authorization: Token ${TOKEN}" \
     --header "Content-Type: text/plain; charset=utf-8" \
     --header "Accept: application/json" \
     --data-binary "$DATA" &> /dev/null

exit $?
