import pandas as pd

print("Scanning....")
df_orders = pd.read_excel('orders-2022-03-21-18-02-34.xlsx')
df_leden = pd.read_excel('Up to date ledenlijst 146 20211116.xlsx')


df_filter = df_orders.drop(df_orders[df_orders.betaalmethode == 'iDEAL'].index)
df_orders2 = df_filter['First Name (Billing)'] + ' ' + df_filter['Last Name (Billing)']
df_leden2 = df_leden['roepnaam'] + ' ' + df_leden.fillna('')['tussenvoegsel'] + ' ' + df_leden['achternaam']
df_leden2 = df_leden2.str.replace('  ', ' ')
df_orders2 = df_orders2.str.normalize('NFKD').str.encode('ascii', errors='ignore').str.decode('utf-8').str.lower()
df_leden2 = df_leden2.str.normalize('NFKD').str.encode('ascii', errors='ignore').str.decode('utf-8').str.lower()


cond = df_orders2.isin(df_leden2)
df_orders2.drop(df_orders2[cond].index, inplace = True)
result = df_orders2.reset_index(drop= True)
print()
print("Mensen die geturfd hebben maar niet in de ledenlijst staan:")
print("-----------------------------------------------------------")
print(result)
