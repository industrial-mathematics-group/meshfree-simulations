#!/bin/bash
# Shows free space and depot usage as seen from a job.
df -h /depot /work /tmp; du -sh /depot/* 2>/dev/null | sort -h; du -sh /depot/artifacts/* 2>/dev/null | sort -h | tail -8
