import requests
from requests.exceptions import HTTPError
import pandas as pd
import sys


def get_protein_information(accessions):
    """
    Get protein information from the UniProt API given a list of accessions.

    Parameters
    ----------
    accessions : list
        List of UniProt accessions as strings.

    Returns
    -------
    pandas.DataFrame
        Dataframe containing the collected protein information.
    """
    print("Getting protein information from UniProt...")
    data = []

    # iterate over accesions, make API requests and store data from returned json object
    for accession in accessions:
        try:
            response = requests.get(f"https://www.ebi.ac.uk/proteins/api/proteins/{accession}")
            response.raise_for_status()

            response_json = response.json()
            protein_info = {
                "Protein Accession": accession,
                "Protein Name": response_json["protein"]["recommendedName"]["fullName"]["value"],
                "Gene": response_json["gene"][0]["name"]["value"],
                "Organism (Scientific)": response_json["organism"]["names"][0]["value"],
                "Organism (Common)": response_json["organism"]["names"][1]["value"],
                "Molecular Weight (Da)": response_json["sequence"]["mass"]
            }
            data.append(protein_info)

        # error handling but no termination of the script since the other accessions might work
        except HTTPError as err:
            print(f"HTTPError for acession {accession}: {err}")
            print("Other accessions will still be processed")
        except Exception as err:
            print(f"Exception for acession {accession}: {err}")
            print("Other accessions will still be processed")

    # create and return result df
    return pd.DataFrame(data)


def get_ensembl_gene_ids(df):
    """
    Get ensembl gene IDs from the Ensembl API and append them to the provided dataframe.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe containing protein information created by get_protein_information().

    Returns
    -------
    pandas.DataFrame
        Modified input dataframe with appended Ensembl gene IDs.
    """
    print("Getting Ensembl gene IDs...")
    ensembl_ids = []

    # iterate over species and genes, make API request and fetch ensembl IDs from returned json object
    for species, gene in zip(df["Organism (Scientific)"], df["Gene"]):
        try:
            response = requests.get(f"https://rest.ensembl.org/xrefs/symbol/{species.lower()}/{gene}",
                                    headers={"Content-Type": "application/json"})
            response.raise_for_status()

            response_json = response.json()
            ensembl_ids.append(response_json[0]["id"])

        # error handling but no termination of the script since the other accessions might work
        except HTTPError as err:
            print(f"HTTPError for species {species} and gene {gene}: {err}")
            print("Other accessions will still be processed")
        except Exception as err:
            print(f"Exception for species {species} and gene {gene}: {err}")
            print("Other accessions will still be processed")

    df["Ensembl Gene ID"] = ensembl_ids
    return df


def get_gene_data(ensembl_ids):
    """
    Get gene information from the Ensembl API given a list of gene IDs.

    Parameters
    ----------
    ensembl_ids : list
        List of Ensembl gene IDs as strings.

    Returns
    -------
    pandas.DataFrame
        Dataframe containing the collected gene information.
    """
    print("Getting gene information from Ensembl...")
    data = []
    headers = {"Content-Type": "application/json", "Accept": "application/json"}

    # make API call and fetch gene data from returned json object
    try:
        response = requests.post("https://rest.ensembl.org/lookup/id", headers=headers, json={"ids": ensembl_ids})
        response.raise_for_status()

        response_json = response.json()
        for ensembl_id in ensembl_ids:
            gene_data = {"Ensembl Gene ID": ensembl_id, "Description": response_json[ensembl_id]["description"],
                         "Seq Region Name": response_json[ensembl_id]["seq_region_name"]}
            data.append(gene_data)

    # error handling resulting in termination since only one API call is made for all genes so no data will be available
    except HTTPError as err:
        print(f"HTTPError when trying to fetch ensembl gene data: {err}")
        print("Terminating program")
        sys.exit()
    except Exception as err:
        print(f"Exception when trying to fetch ensembl gene data: {err}")
        print("Terminating program")
        sys.exit()

    return pd.DataFrame(data)


# run the functions
if __name__ == "__main__":
    accessions = ["P12345", "Q8N726", "O00255"]
    prot_inf_df = get_protein_information(accessions)
    prot_inf_df = get_ensembl_gene_ids(prot_inf_df)
    gene_inf_df = get_gene_data(prot_inf_df["Ensembl Gene ID"].tolist())

    # merge dataframes
    result = prot_inf_df.merge(gene_inf_df, how="outer", on="Ensembl Gene ID")

    # store result df as excel file
    outfile = "protein_gene_analysis.xlsx"
    print("Writing results to " + outfile)
    result.to_excel(outfile, index=False)
