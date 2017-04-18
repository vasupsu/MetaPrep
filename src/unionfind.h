uint32_t findRoot1 (uint32_t *parent, uint32_t elem)
{
	uint32_t root=elem;
	while (parent[root] != root)
	{
		root = parent[root];
//		*numLoads = *numLoads+1;
	}
	while (elem != root)
	{
		uint32_t curr_node = parent[elem];
		parent[elem] = root;
		elem = curr_node;
//		*numLoads = *numLoads+1;
//		*numStores = *numStores+1;
	}
	return elem;
}
uint32_t findRoot (uint32_t *parent, uint32_t elem, int *numRandAccess)
{
        while (parent[elem] != elem)
        {
		//onepass
                parent[elem] = parent[parent[elem]];
                elem = parent[elem];
                *numRandAccess = *numRandAccess + 1;
//		*numLoads = *numLoads+2;
        }
        return elem;
}

//Assume u<v
void unionOp (uint32_t *parent, uint32_t *sizes, uint32_t u, uint32_t v)
{
	//simple
	parent[v]=u;
//	fprintf (fp, "UF %u(size %lu)->%u (size %lu)\n", v, sizes[v], u, sizes[u]);
//	sizes[u]+=sizes[v];

	//weighted
/*	if (sizes[u] < sizes[v])
	{
		parent[u]=v;
		sizes[v]+=sizes[u];
	}
	else
	{
		parent[v]=u;
		sizes[u]+=sizes[v];
	}*/
}
