#!/bin/bash
until [ -s "$1" ] && grep -q -E "$2" "$1"; do sleep 20; done
tail -c 600 "$1"
