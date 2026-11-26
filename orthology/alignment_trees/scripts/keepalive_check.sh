#!/bin/bash
ps -eo pid,etime,cmd | grep -E "sleep 600|sleep infinity" | grep -v grep
uptime
