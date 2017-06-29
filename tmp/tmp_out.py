channel_out = open('/Users/mojarro/Desktop/carrierseq_channels_clean.lst', 'r')
channel_list = channel_out.readlines()

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

# Next we create a dictionary where each element is in the format of "element: frequency"
channel_freq = {x:channel_list_num.count(x) for x in channel_list_num}

print channel_freq