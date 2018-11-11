if test_detection_probs is None:
    test_detection_probs = [[0.0] * 20 + [0.8] * 40 + [0.2] * 10 + [0.0] * 20]
if false_pos_prob is None:
    false_pos_prob = [0] * len(test_detection_probs)

store_test_results = [[None] * len(days_since_infection) for n in range(len(test_detection_probs))]

for tt in range(len(test_detection_probs)):
    cur_test_probs = test_detection_probs[tt]
    for ii in range(len(days_since_infection)):
        if days_since_infection[ii] >= len(cur_test_probs):
            cur_prob = 0
        else:
            cur_prob = cur_test_probs[days_since_infection[ii]]
        # could test positive from actual infection or from false positive
        prob_pos = 1 - (1 - cur_prob) * (1 - false_pos_prob[tt])
        store_test_results[tt][ii] = np.random.binomial(n=1, p=prob_pos)