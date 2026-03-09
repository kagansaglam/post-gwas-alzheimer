def check_stock(inventory_list,product_name):
    inventory_list=inventory
    product_name=int(input('write a name:'))
    for i in product_name:
        if i == inventory_list:
            print(i)
        else:
            print('product not found')
check_stock(inventory_list=inventory,product_name='laptop')