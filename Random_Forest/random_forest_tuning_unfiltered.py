import pandas as pd
import os
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import r2_score
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Load your data into a pandas DataFrame
df = pd.read_csv('out_subelements.csv')
subset = pd.read_csv(os.path.join("subset.csv"))

# Merge dataframes based on "genome" column
merged_data = pd.merge(df, subset, on="genome", how="left")

# Replace NAs in "rank" column of df with values from subset
merged_data.loc[~merged_data['rank_y'].isna(), 'rank_x'] = merged_data.loc[~merged_data['rank_y'].isna(), 'rank_y']

# Remove unnecessary columns
merged_data = merged_data.drop(columns=['rank_y'])

# Rename columns if needed
merged_data = merged_data.rename(columns={'rank_x': 'rank'})

# Filter columns based on Stats_filtered criteria
Stats_filtered = pd.read_table(os.path.join("out_subelements.key.txt"), sep="\t")


# Filter columns in df based on common subelements
common_columns = list(set(Stats_filtered['subelement']) & set(df.columns))
columns_to_keep = ['genome', 'rank'] + common_columns
df = df[columns_to_keep]

# Split the data into features (X) and target variable (y)
X = df.drop(['genome', 'rank'], axis=1)  # Features
y = df['rank']  # Target variable

# Define hyperparameter grid for tuning
param_grid = {
    'n_estimators': [50, 100, 200],        # Number of trees
    'max_depth': [None, 10, 20, 30],      # Maximum depth of trees
    'min_samples_split': [2, 5, 10],      # Minimum samples required to split an internal node
    'min_samples_leaf': [1, 2, 4]         # Minimum samples required to be at a leaf node
}

# Initialize an empty DataFrame to store the results
results_data = []

# Perform hyperparameter tuning with GridSearchCV
grid_search = GridSearchCV(estimator=RandomForestRegressor(random_state=999),
                           param_grid=param_grid,
                           scoring='r2',
                           cv=5,  # 5-fold cross-validation
                           verbose=2,
                           n_jobs=-1)

# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=50, train_size=207, random_state=999)

# Fit the model using GridSearchCV
grid_search.fit(X_train, y_train)

# Get the best parameters from GridSearchCV
best_params = grid_search.best_params_
print("Best Hyperparameters:", best_params)

# Train the best model on the entire training set
best_rf_regressor = RandomForestRegressor(**best_params, random_state=999, n_jobs= -1)
best_rf_regressor.fit(X_train, y_train)

# Predict ranks for the training set
y_train_pred = best_rf_regressor.predict(X_train)

# Predict ranks for the test set
y_test_pred = best_rf_regressor.predict(X_test)

# Calculate R-squared for the training set
r2_train = r2_score(y_train, y_train_pred)

# Calculate R-squared for the test set
r2_test = r2_score(y_test, y_test_pred)

# Append results to DataFrame
results_data.append({'Random State': 999, 'Train R-squared': r2_train, 'Test R-squared': r2_test})

# Create DataFrame from the list of results
results = pd.DataFrame(results_data)

# Save results to CSV
results.to_csv('random_forest_results_unfiltered.csv', index=False)


# Initialize the Random Forest regressor with the best hyperparameters
#best_params = {'max_depth': 10, 'min_samples_leaf': 1, 'min_samples_split': 10, 'n_estimators': 50}
best_rf_regressor = RandomForestRegressor(**best_params, random_state=999)

# Train the model on the entire dataset
best_rf_regressor.fit(X, y)

# Get feature importances
feature_importances = best_rf_regressor.feature_importances_

# Create a DataFrame for feature importances
feature_importance_df = pd.DataFrame({
    'Feature': X.columns,
    'Importance': feature_importances
})

# Sort the DataFrame by absolute value of importance
feature_importance_df['AbsImportance'] = feature_importance_df['Importance'].abs()
feature_importance_df = feature_importance_df.sort_values(by='AbsImportance', ascending=False)

# Display the top 30 most important features by absolute value of importance
top_30_features = feature_importance_df.head(30)
print(top_30_features)

# Save the top 30 feature importances to a CSV file
top_30_features.to_csv('top_30_feature_importances.csv', index=False)

# Count the number of features with non-zero importance
num_important_features = np.sum(feature_importance_df['Importance'] > 0)
print(f"Number of features with non-zero importance: {num_important_features}")
print(f"Number of strains in the training dataset: {X_train.shape[0]}")
print(f"Number of strains in the test dataset: {X_test.shape[0]}")
# Assuming you have the following from your model
# y_test: actual values for the test set
# y_test_pred: predicted values from your model

# Create a figure with two subplots (1 row, 2 columns)
fig, axes = plt.subplots(1, 2, figsize=(15, 6))  # Adjust figsize as needed

# Create the first scatter plot (Figure 1A)
sns.scatterplot(x=y_test, y=y_test_pred, ax=axes[1])
axes[1].plot([0, 1], [0, 1], 'r--', lw=2)  # 45-degree line from (0,0) to (1,1)
axes[1].set_xlabel('Actual Values')
axes[1].set_ylabel('Predicted Values')
axes[1].set_title('Figure 1B: Predicted vs Actual Values (Test)')
axes[1].set_xlim(0.0, 1.0)  # Set x-axis limit
axes[1].set_ylim(0.0, 1.0)  # Set y-axis limit
axes[1].grid()

# Add R² value to the first plot
axes[1].text(0.05, 0.95, f'R² = {r2_test:.2f}', fontsize=12, ha='left', va='top', transform=axes[1].transAxes)

# Create the second scatter plot (Figure 1B)
sns.scatterplot(x=y_train, y=y_train_pred, ax=axes[0])
axes[0].plot([0, 1], [0, 1], 'r--', lw=2)  # 45-degree line from (0,0) to (1,1)
axes[0].set_xlabel('Actual Values')
axes[0].set_ylabel('Predicted Values')
axes[0].set_title('Figure 1A: Predicted vs Actual Values (Train)')
axes[0].set_xlim(0.0, 1.0)  # Set x-axis limit
axes[0].set_ylim(0.0, 1.0)  # Set y-axis limit
axes[0].grid()

# Add R² value to the second plot
axes[0].text(0.05, 0.95, f'R² = {r2_train:.2f}', fontsize=12, ha='left', va='top', transform=axes[0].transAxes)

# Save the figure with both subplots
plt.savefig('predicted_vs_actual_filtered_combined.png', dpi=300, bbox_inches='tight')
plt.show()

# Show the plot
