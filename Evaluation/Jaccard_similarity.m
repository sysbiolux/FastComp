function sim = Jaccard_similarity(origin,pred)

sim = ((sum(origin ~= pred) + sum(origin == pred)) - sum(origin == pred)) / (sum(origin ~= pred) + sum(origin == pred));

    