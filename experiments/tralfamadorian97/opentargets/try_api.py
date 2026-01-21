import requests
import json

# NOTE: The endpoint is api.GENETICS, not api.platform
base_url = "https://api.genetics.opentargets.org/graphql"

# The 'manhattan' query is a root-level field in the Genetics schema
query = """
query WindowQuery($studyId: String!, $pageIndex: Int!, $pageSize: Int!) {
  manhattan(studyId: $studyId, pageIndex: $pageIndex, pageSize: $pageSize) {
    associations {
      variant {
        id
        rsId
        position
        chromosome
      }
      beta
      se
      pval
    }
  }
}
"""

# Example ID for a UKBB PPP protein (IL6)
variables = {
    "studyId": "GCST90266624",
    "pageIndex": 0,
    "pageSize": 100
}

response = requests.post(base_url, json={"query": query, "variables": variables})
data = response.json()

print(json.dumps(data, indent=2))