from HorusAPI import PluginBlock, PluginVariable, VariableTypes

from pandadock.protein import Protein
from pandadock.ligand import Ligand
from pandadock.scoring import CompositeScoringFunction, EnhancedScoringFunction
from pandadock.search import GeneticAlgorithm
from pandadock.utils import save_docking_results
from pandadock.preparation import prepare_protein, prepare_ligand


# Inputs of the block
input_protein = PluginVariable(
    id="input_protein",
    name="Input protein",
    description="The target protein (PDB or MOL2)",
    type=VariableTypes.FILE,
    allowedValues=["pdb", "mol2"],
)

input_ligand = PluginVariable(
    id="input_ligand",
    name="Ligand",
    description="Ligand to dock into the target protein",
    type=VariableTypes.FILE,
    allowedValues=["pdb", "sdf", "mol"],
)

sphere_center = PluginVariable(
    id="dock_center",
    name="Docking region",
    description="A sphere denoting the dock search region",
    type=VariableTypes.SPHERE,
)

# Block internal variables
search_algorithm_scoring_function = PluginVariable(
    id="search_algorithm_scoring_function",
    name="Scoring Function",
    description="The scoring function used for docking",
    type=VariableTypes.STRING_LIST,
    allowedValues=["Enhanced", "Composite"],
    defaultValue="Enhanced",
    category="Docking options",
)

search_algorithm_max_iterations = PluginVariable(
    id="search_algorithm_max_iterations",
    name="Max Iterations",
    description="The maximum number of iterations for the genetic algorithm",
    type=VariableTypes.NUMBER,
    defaultValue=1000,
    category="Docking options",
)

search_algorithm_population_size = PluginVariable(
    id="search_algorithm_population_size",
    name="Population Size",
    description="The size of the population in the genetic algorithm",
    type=VariableTypes.NUMBER,
    defaultValue=100,
    category="Docking options",
)

search_algorithm_mutation_rate = PluginVariable(
    id="search_algorithm_mutation_rate",
    name="Mutation Rate",
    description="The mutation rate for the genetic algorithm",
    type=VariableTypes.NUMBER,
    defaultValue=0.3,
    category="Docking options",
)


add_hydrogens_protein = PluginVariable(
    id="add_hydrogens_protein",
    name="Add Hydrogens",
    description="Wether to protonate or not the protein",
    type=VariableTypes.BOOLEAN,
    defaultValue=True,
    category="Protein",
)
ph_protein = PluginVariable(
    id="ph_protein",
    name="pH",
    description="pH used to protonate the protein. Only used when 'Add Hydrogens' is set to true",
    type=VariableTypes.NUMBER_RANGE,
    allowedValues=[0.1, 14, 0.1],
    defaultValue=7.4,
    category="Protein",
)

add_hydrogens_ligand = PluginVariable(
    id="add_hydrogens_ligand",
    name="Add Hydrogens",
    description="Wether to protonate or not the ligand",
    type=VariableTypes.BOOLEAN,
    defaultValue=True,
    category="Ligand",
)

# Outputs of the block
docking_results = PluginVariable(
    id="docking_results",
    name="Docking output",
    description="The folder containing the docked results",
    type=VariableTypes.FOLDER,
)


def pandadock(block: PluginBlock):

    protein_input_path = block.inputs[input_protein.id]
    ligand_input_path = block.inputs[input_ligand.id]
    docking_sphere_value = block.inputs[sphere_center.id]

    docking_center = [p for p in docking_sphere_value["center"].values()]
    docking_radius = docking_sphere_value["radius"]

    add_hydrogens_protein_value = block.variables[add_hydrogens_protein.id]
    ph_protein_value = block.variables[ph_protein.id]
    add_hydrogens_ligand_value = block.variables[add_hydrogens_ligand.id]

    search_algorithm_scoring_function_value = block.variables[
        search_algorithm_scoring_function.id
    ]
    search_algorithm_max_iterations_value = block.variables[
        search_algorithm_max_iterations.id
    ]
    search_algorithm_population_size_value = block.variables[
        search_algorithm_population_size.id
    ]
    search_algorithm_mutation_rate_value = block.variables[
        search_algorithm_mutation_rate.id
    ]

    # Prepare molecules (optional)
    prepared_protein = prepare_protein(
        protein_input_path,
        add_hydrogens=add_hydrogens_protein_value,
        ph=ph_protein_value,
    )
    prepared_ligand = prepare_ligand(
        ligand_input_path,
        minimize=True,
        add_hydrogens=add_hydrogens_ligand_value,
    )

    # Load protein and ligand
    protein = Protein(prepared_protein)
    ligand = Ligand(prepared_ligand)

    # Define active site (optional)
    protein.define_active_site(docking_center, docking_radius)

    # Create scoring function (choose basic or enhanced)
    scoring_function = None  # or CompositeScoringFunction()

    if search_algorithm_scoring_function_value == "Composite":
        scoring_function = CompositeScoringFunction()
    elif search_algorithm_scoring_function_value == "Enhanced":
        scoring_function = EnhancedScoringFunction()
    else:
        raise ValueError(
            f"Unsupported scorting function '{search_algorithm_scoring_function_value}'"
        )

    # Create search algorithm
    search_algorithm = GeneticAlgorithm(
        scoring_function,
        max_iterations=search_algorithm_max_iterations_value,
        population_size=search_algorithm_population_size_value,
        mutation_rate=search_algorithm_mutation_rate_value,
    )

    # Perform docking
    results = search_algorithm.search(protein, ligand)

    # Apply local optimization (optional)
    optimized_results = []
    for i, (pose, score) in enumerate(sorted(results, key=lambda x: x[1])[:5]):
        opt_pose, opt_score = search_algorithm._local_optimization(pose, protein)
        optimized_results.append((opt_pose, opt_score))

    # Save results
    save_docking_results(optimized_results, "docking_results")


pandadock_block = PluginBlock(
    id="pandadock",
    name="PandaDock",
    description="Perform a PandaDock docking",
    inputs=[input_protein, input_ligand, sphere_center],
    variables=[
        search_algorithm_scoring_function,
        search_algorithm_population_size,
        search_algorithm_max_iterations,
        search_algorithm_mutation_rate,
        add_hydrogens_protein,
        ph_protein,
        add_hydrogens_ligand,
    ],
    outputs=[docking_results],
    action=pandadock,
)
