
using NetworkInference

nodes_healthy = get_nodes("./CD4network/Healthy.txt")
nodes_mild = get_nodes("./CD4network/Mild.txt")
nodes_severe = get_nodes("./CD4network/Severe.txt")

inferred_network_healthy = InferredNetwork(PIDCNetworkInference(), nodes_healthy)
inferred_network_mild = InferredNetwork(PIDCNetworkInference(), nodes_mild)
inferred_network_severe = InferredNetwork(PIDCNetworkInference(), nodes_severe)

write_network_file("./CD4network/HealthyPIDC.txt", inferred_network_healthy)
write_network_file("./CD4network/MildPIDC.txt", inferred_network_mild)
write_network_file("./CD4network/SeverePIDC.txt", inferred_network_severe)


Threads.@threads for i in 1:500

	nodes_HM_healthy = get_nodes("./CD4network/HealthyMildPermDat/permH$i.txt")
	inferred_network_HM_healthy = InferredNetwork(PIDCNetworkInference(), nodes_HM_healthy)
	write_network_file("./CD4network/HealthyMildPIDC/netH$i.txt", inferred_network_HM_healthy)
	
	nodes_HM_mild = get_nodes("./CD4network/HealthyMildPermDat/permM$i.txt")
	inferred_network_HM_mild = InferredNetwork(PIDCNetworkInference(), nodes_HM_mild)
	write_network_file("./CD4network/HealthyMildPIDC/netM$i.txt", inferred_network_HM_mild)
	
	
	nodes_HS_healthy = get_nodes("./CD4network/HealthySeverePermDat/permH$i.txt")
	inferred_network_HS_healthy = InferredNetwork(PIDCNetworkInference(), nodes_HS_healthy)
	write_network_file("./CD4network/HealthySeverePIDC/netH$i.txt", inferred_network_HS_healthy)
	
	nodes_HS_severe = get_nodes("./CD4network/HealthySeverePermDat/permS$i.txt")
	inferred_network_HS_severe = InferredNetwork(PIDCNetworkInference(), nodes_HS_severe)
	write_network_file("./CD4network/HealthySeverePIDC/netS$i.txt", inferred_network_HS_severe)
	
	
	nodes_MS_mild = get_nodes("./CD4network/MildSeverePermDat/permM$i.txt")
	inferred_network_MS_mild = InferredNetwork(PIDCNetworkInference(), nodes_MS_mild)
	write_network_file("./CD4network/MildSeverePIDC/netM$i.txt", inferred_network_MS_mild)
	
	nodes_MS_severe = get_nodes("./CD4network/MildSeverePermDat/permS$i.txt")
	inferred_network_MS_severe = InferredNetwork(PIDCNetworkInference(), nodes_MS_severe)
	write_network_file("./CD4network/MildSeverePIDC/netS$i.txt", inferred_network_MS_severe)
	
	
	println(i)
	
end