input = 'A*02:01:01#A*02:06:01'
input_array = input.split('#')
h_type =[]
for i in range(len(input_array)):
    start = input_array[i].index('*')
    end = input_array[i].rindex(':')
    type = input_array[i][start:end]
    h_type.append(type)
print('/'.join(h_type))