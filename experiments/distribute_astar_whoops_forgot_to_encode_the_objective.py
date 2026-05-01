"""
Here, we will implement a distributed A* algorithm that can handle
multiple configurations for a breeding program optimization problem. The algorithm will be designed to run in parallel across different instances, and we will visualize the results using box plots to compare the performance of various configurations.
"""

import time
import random
import eugene_rs
import eugene.utils as eu
import pymongo as pm

LOCI = list(range(2, 17))

INSTANCES = {
    n_loci: [eu.random_distribute_instance(n_loci) for _ in range(100)]
    for n_loci in LOCI
}

# Connect to MongoDB
client = pm.MongoClient("mongodb://localhost:27017/")

# Create or get the database
db = client["eugene_breeding_programs"]

# Uncomment the line below to clean up the collection if it exists
# db.drop_collection("distribute_astar_results")

# Create or get the collection
collection = db["distribute_astar_results"]


# print the contents of the collection
def print_collection_contents():
    for document in collection.find():
        print(document)


# Uncomment to print the contents of the collection
# print_collection_contents()

for document in collection.find():
    print(document)
    if "objective" in document:
        print(f"Skipping document with _id={document['_id']} due to existing objective")
        continue
    output = (
        eugene_rs.min_cross.distribute_astar.breeding_program_distribute_general_python(
            document["instance"],
            300,
            diving=document["diving"],
            full_join=document["full_join"],
            dominance=document["dominance"],
        )
    )
    print(f"Output: {output}")
    # Store the results in the collection
    if output is None:
        collection.update_one(
            {"_id": document["_id"]},
            {"$set": {"objective": None, "children_created_by_branching": None}},
        )
    else:
        collection.update_one(
            {"_id": document["_id"]},
            {
                "$set": {
                    "objective": output["objective"],
                    "children_created_by_branching": output[
                        "children_created_by_branching"
                    ],
                }
            },
        )
