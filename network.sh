#!/bin/bash

# Network interface name, replace with your actual interface name
INTERFACE="lo"

# Bandwidth limit
BANDWIDTH="100mbit"

# Delay
DELAY="40ms"

# Setup the root qdisc
sudo tc qdisc add dev $INTERFACE root handle 1: htb default 10

# Add class for rate limiting
sudo tc class add dev $INTERFACE parent 1: classid 1:1 htb rate $BANDWIDTH

# Add qdisc for delay
sudo tc qdisc add dev $INTERFACE parent 1:1 handle 10: netem delay $DELAY

echo "Network settings applied: $BANDWIDTH bandwidth limit and $DELAY delay on $INTERFACE."
