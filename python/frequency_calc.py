import sys

channel_out = open(sys.argv[1], 'r')
channel_list = channel_out.readlines()

xcrit_value_txt = open(sys.argv[2], 'r')
xcrit_value = xcrit_value_txt.read()
xcrit = float(xcrit_value)

# channel_list: A list containing strings of each channel with newline characters
# new_channel_list: A list containing strings of each channel
# channel_list_num: A list containing integers of each channel
 
# First we strip the newline characters in a loop
new_channel_list = []
i = 0
for element in channel_list:
 	new_channel_list.append(channel_list[i].rstrip('\n'))
 	i = i + 1
 
# Next we convert each element to an integer in a loop
ind = 0
channel_list_num = []
for element in new_channel_list:
 	channel_list_num.append(int(new_channel_list[ind]))
 	ind = ind + 1
 
# Next we create a dictionary where each element is in the format of "channel: frequency"
channel_freq = {x:channel_list_num.count(x) for x in channel_list_num}

# print channel_freq 
target_channels = dict()
hqnr_channels = dict()

for channel in channel_freq:
	if channel_freq[channel] <= xcrit:
		target_channels[channel] = channel_freq[channel]
	else:
		hqnr_channels[channel] = channel_freq[channel]

# Save roi channel frequency dictionary
with open(sys.argv[3], 'w') as f:
    sys.stdout = f
    print channel_freq
    
# save hqnr channel dictionary
with open(sys.argv[4], 'w') as f:
    sys.stdout = f
    print hqnr_channels
    
# save target reads channel dictionary
with open(sys.argv[5], 'w') as f:
    sys.stdout = f
    print target_channels
 
# print only target channels used for sorting    
with open(sys.argv[6], 'w') as f:
    sys.stdout = f
    for item in target_channels.keys():
    	print item
    