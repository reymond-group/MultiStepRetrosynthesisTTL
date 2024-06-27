

from tqdm import tqdm
import pandas as pd
import xml


from Reaxys_API import Reaxys_API
api = Reaxys_API()
api.connect(
    url="https://www.reaxys.com/reaxys/api",
    url_main='',
	username="yourusername",
	password="yourpassword",
	callername="yourcallername")

if api.sessionid == '': 
    print('Unable to connect, exiting..')
    exit()
print('Connected, session ID:', api.sessionid)

list_IDs = []
with open("ENZR_Rxn_IDs.txt", "r") as f:
    for line in f:
        list_IDs.append(line)
        

field_names = ["Reaction_ID", "Reactant_XRN", "Reactant", "Product_XRN", "Product", "Number_of_Reaction_Details", "Registration_String", "Reaction_Rank", 
"Maximum_Yield", "Reaction_Structure_Keywords", "Record_Type", "Reactant_Availability", "Product_Availability", "Maximum_Publication_Year", 
"Number_of_References", "Maximum_Product_Molweight", "Reaction_Entry_Date", "Reaction_Update_Date", "Citation_Pointer", "Reference_(Citation)", 
"Source", "Reaction_Classification", "Reaction_Score", "Number_of_Reaction_Steps", "Multi-step_Scheme", "Multi-step_Details", 
"Yield", "Yield_(numerical)", "Yield_(optical)", "Reagent_XRN", "Reagent", "Catalyst_XRN", "Catalyst", "Solvent_XRN", "Solvent", 
"Reagent/Catalyst", "Other_Conditions", "Entry_Date", "Reaction_mol", "Reactant_mol", "Product_mol"]

fields = ["RX.ID", "RX.RXRN", "RX.RCT", "RX.PXRN", "RX.PRO", "RX.NVAR", "RX.REG", "RX.RANK", "RX.MYD", "RX.SKW", "RX.RTYP", "RX.RAVAIL", 
"RX.PAVAIL", "RX.MAXPUB", "RX.NUMREF", "RX.MAXPMW", "RX.ED", "RX.UPD", "RXD.L", "RXD.citation", "RXD.SOURCE", "RXD.CL", "RXD.SCO", 
"RXD.STP", "RXD.MID", "RXD.MTEXT", "RXD.YD", "RXD.NYD", "RXD.YDO", "RXD.RGTXRN", "RXD.RGT", "RXD.CATXRN", "RXD.CAT", "RXD.SOLXRN", 
"RXD.SOL", "RXD.RGTCAT", "RXD.COND", "RXD.DED", "RY.STR", "RY.RCT", "RY.PRO"]

REFS = {
    "RX.RXRN":      'RX01', 
    "RX.RCT":       'RX01', 
    "RX.PXRN":      'RX02', 
    "RX.PRO":       'RX02', 
    "RXD.RGTXRN":   'RXD03', 
    "RXD.RGT":      'RXD03', 
    "RXD.CATXRN":   'RXD04', 
    "RXD.CAT":      'RXD04', 
    "RXD.SOLXRN":   'RXD05', 
    "RXD.SOL":      'RXD05'
    }

FIELDS_FORCE_REFERENCE_SEARCH = ["RX.RXRN", "RX.RCT", "RX.PXRN", "RX.PRO", "RXD.RGTXRN", "RXD.RGT", "RXD.CATXRN", "RXD.CAT", "RXD.SOLXRN", "RXD.SOL"]

def __cleaner_get_field(element, highlight_only=False):
    # Concatenate text values if highlight is present
    if element.getAttribute('highlight') == 'true':
        return ''.join([e.data
                        if type(e) == xml.dom.minidom.Text
                        else e.childNodes[0].data for e in element.childNodes])

    # If node contains further sub-nodes: return full xml.
    elif len(element.childNodes) > 1 and highlight_only is False:
        try:
            return ''.join([e.data
                            if type(e) == xml.dom.minidom.Text
                            else e.childNodes[0].data for e in element.childNodes])
        except:
            return element.toxml()

    # If node does not conatin further sub-nodes: return node value.
    elif len(element.childNodes) == 1 and highlight_only is False:
        return element.childNodes[0].nodeValue
    
    return ''

def get_field_content_rxn_fields_corrected(response_xml, fields, highlight_only=False):
    
    response_dom = xml.dom.minidom.parseString(response_xml)

    all_responses_rxn = []
    # For each reaction in the query batch:
    for response_dom_rxn in response_dom.getElementsByTagName('reaction'):
        
        rxn_respond = []
        # For each field we want to extract, depending on the requirements:
        for field in fields:
            # If part of exceptions, then the extraction proceed per corresponding reference and reports a list. 
            if field in FIELDS_FORCE_REFERENCE_SEARCH:
                
                # Then, should extract each element of the corresponding reference
                reference = REFS[field]
                
                ref_list = []
                for ref in response_dom_rxn.getElementsByTagName(reference):
                    if len(ref.getElementsByTagName(field)) == 1:       ref_list.append(__cleaner_get_field(ref.getElementsByTagName(field)[0], highlight_only=highlight_only))
                    elif len(ref.getElementsByTagName(field)) == 0:     ref_list.append('')
                    else:
                        print('ERROR: More than one field found for the reference: ', reference, ' and field: ', field)
                        ref_list.append('')
                
                rxn_respond.append(ref_list)
                    
            else:
                field_content = []
                for element in response_dom_rxn.getElementsByTagName(field):
                    field_content.append(__cleaner_get_field(element, highlight_only=highlight_only))
                    
                rxn_respond.append(field_content)
        
        all_responses_rxn.append(rxn_respond)
    return all_responses_rxn


MIN = 0
MAX = len(list_IDs)
STEP = 1000 # Cannot exceed 1000 IDs per request

df = pd.DataFrame(columns=field_names)

for i in tqdm(range(MIN, MAX, STEP)):
    where_clause = "RX.ID = " + ' OR RX.ID = '.join(list_IDs[i:i+STEP])
    api.select(
        dbname="RX", 
        context="R", 
        where_clause=where_clause, 
        order_by="",
        options="WORKER,NO_CORESULT"
        )
    
    ide_yy_response = api.retrieve(
        resultname=api.resultname, 
        select_items=['RX', 'RXD', 'RY'],
        first_item=1,   # No 0 in Reaxys
        last_item=STEP, 
        order_by="", 
        group_by="", 
        group_item="",
        options="HITONLY,EXPORT=true,ISSUE_RXN=true,COMPRESS=false,OMIT_MAPS=true", 
        dbname='RX', context=None
        )
    
    reply = get_field_content_rxn_fields_corrected(
        response_xml=ide_yy_response,
        fields=fields,
        highlight_only=False, 
        )
    
    for rxn in range(0, len(reply)):
        for field in range(0, len(reply[rxn])):
            df.at[i+rxn, field_names[field]] = reply[rxn][field]
    

for column in tqdm(df.columns):
    column_max_list_len = 0
    
    # Check the max length of the list in the column:
    for el in range(0, len(df)):
        if isinstance(df.at[el, column], list):
            if len(df.at[el, column]) > column_max_list_len:
                column_max_list_len = len(df.at[el, column])
        else:
            continue
    
    if column_max_list_len == 1:
        for el in range(0, len(df)):
            if len(df.at[el, column]) == 1:
                df.at[el, column] = df.at[el, column][0]
            else:
                df.at[el, column] = pd.NA
                
                
df.to_pickle('ENZR_Dataset_raw.pkl')