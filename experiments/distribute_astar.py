import time
import random
import eugene_rs
import eugene.utils as eu
import pymongo as pm

random.seed(0)

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
print_collection_contents()

CONFIGS = [
    {"diving": True, "full_join": True, "dominance": False},
    {"diving": False, "full_join": True, "dominance": False},
    {"diving": False, "full_join": False, "dominance": False},
    {"diving": False, "full_join": True, "dominance": True},
    {"diving": False, "full_join": False, "dominance": True},
]

TITLES = [
    "Distribute A* with full-joins and diving:",
    "Distribute A* with full_joins:",
    "Distribute A* no full-joins or diving",
    "Distribute A* with dominance and full-joins:",
    "Distribute A* with dominance but no full-joins or diving",
]

for conf, title in zip(CONFIGS, TITLES):
    print(title)
    for n_loci in LOCI:
        for instance_number, instance in enumerate(INSTANCES[n_loci]):
            # Check if the instance already exists in the collection
            if collection.find_one(
                {"n_loci": n_loci, "instance_number": instance_number, **conf}
            ):
                print(
                    f"Skipping n_loci={n_loci}, instance={instance} with config={conf} (already exists)"
                )
                continue

            # If the instance does not exist, run the algorithm and store the results
            print(
                f"Running for n_loci={n_loci}, instance={instance}, config={conf}"
            )
            start = time.time()
            for _ in range(8):
                output = eugene_rs.min_cross.distribute_astar.breeding_program_distribute_general_python(
                    instance, 300, **conf
                )
            time_measured = (time.time() - start) / 8
            if output is None:
                print(
                    f"Output is None for n_loci={n_loci}, instance={instance}"
                )
                collection.insert_one(
                    {
                        "n_loci": n_loci,
                        "instance": instance,
                        "instance_number": instance_number,
                        "time": time_measured,
                        "expansions": None,
                        "pushed_nodes": None,
                        "objective": None,
                        "timeout": True,
                        **conf,
                    }
                )
            else:
                collection.insert_one(
                    {
                        "n_loci": n_loci,
                        "instance": instance,
                        "instance_number": instance_number,
                        "time": time_measured,
                        "expansions": output["expansions"],
                        "pushed_nodes": output["pushed_nodes"],
                        "objective": output["objective"],
                        "timeout": False,
                        **conf,
                    }
                )
