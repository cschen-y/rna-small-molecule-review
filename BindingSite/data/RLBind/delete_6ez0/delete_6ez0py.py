import pickle


pdb_name = '6ez0'

def delete_node():
    with open('../test18_all.pkl', 'rb') as f:
        data = pickle.load(f)
    with open('../mot11_T18.pkl', 'rb') as f:
        mot11_data = pickle.load(f)

    length = 0
    all_node = []
    all_mot11 = []
    flat = False
    for node,mot11 in zip(data, mot11_data):
        if node[3] != pdb_name and flat is False:
            all_node.append(node)
            all_mot11.append(mot11)
        elif flat is True and node[3] != pdb_name:
            node_list = list(node)
            node_list[0] = node_list[0] - length
            node_list[2] = node_list[2] - 1
            all_node.append(tuple(node_list))
            all_mot11.append(mot11)
        else:
            length = int(node[4])
            flat = True

    with open('./delete_data/test18_all.pkl', 'wb') as f:
        pickle.dump(all_node , f)
    with open('./delete_data/mot11_T18.pkl', 'wb') as f:
        pickle.dump(all_mot11 , f)
    with open('./delete_data/test18_index.pkl', 'wb') as f:
        train_index = list(range(len(all_node)))
        pickle.dump(train_index , f)
    

def delete_data():
    data_all = []
    label_all = []
    with open('../data_T18.pkl', 'rb') as f:
        data = pickle.load(f)
    with open('../label_T18.pkl', 'rb') as f:
        label = pickle.load(f)
    for i in range(len(data)):
        print(i)
        if i == 10:
            continue
        else:
            label_all.append(label[i])
            data_all.append(data[i])
    with open('./delete_data/data_T18.pkl', 'wb') as f:
        pickle.dump(data_all , f)
    with open('./delete_data/label_T18.pkl', 'wb') as f:
        pickle.dump(label_all , f)
    


delete_node()
delete_data()


