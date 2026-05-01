import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sn
from pymongo import MongoClient


# Connect to MongoDB
client = MongoClient("mongodb://localhost:27017/")

# Create or get the database
db = client["eugene_breeding_programs"]

# Create or get the collection
collection = db["distribute_astar_results"]

# Create a DataFrame from the MongoDB collection
df = pd.DataFrame(list(collection.find()))
print(df.info())

exit()

CONFIGS = [
    {"diving": True, "full_join": True, "dominance": False},
    {"diving": False, "full_join": True, "dominance": False},
    {"diving": False, "full_join": False, "dominance": False},
    {"diving": False, "full_join": True, "dominance": True},
    {"diving": False, "full_join": False, "dominance": True},
]

TITLES = [
    "Distribute A* with full-joins and diving",
    "Distribute A* with full_joins",
    "Distribute A* no full-joins or diving",
    "Distribute A* with dominance and full-joins",
    "Distribute A* with dominance but no full-joins or diving",
]


def plot_results(df, config, title):
    filtered_df = df[(df['diving'] == config['diving']) &
                     (df['full_join'] == config['full_join']) &
                     (df['dominance'] == config['dominance'])]

    plt.figure(figsize=(10, 6))
    sn.boxplot(x='n_loci', y='time', data=filtered_df)
    plt.title(title)
    plt.xlabel('Number of Loci')
    plt.ylabel('Time (seconds)')
    plt.xticks(rotation=45)
    plt.yscale('log')  # Use logarithmic scale for better visibility
    plt.tight_layout()
    plt.show()


for conf, title in zip(CONFIGS, TITLES):
    print(title)
    plot_results(df, conf, title)


# Plot the distribution of objectives for all configurations and varying n_loci
def plot_objective_distribution(df):
    plt.figure(figsize=(12, 8))
    sn.boxplot(x='n_loci', y='objective', hue='config', data=df.loc[(df["full_join"] == False) & (df["dominance"] == False)], palette='Set2')
    plt.title('Distribution of Objectives by Number of Loci and Configuration')
    plt.xlabel('Number of Loci')
    plt.ylabel('Objective Value')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()


plot_objective_distribution(df)
