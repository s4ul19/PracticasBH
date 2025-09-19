function andpocalipse(cells)
  local tab = {}
  for i,v in pairs(cells) do
    tab[v] = true
  end
  for i=1,15 do
    if tab[i] then
      tex.sprint("\\cellcolor[HTML]{C0C0C0}&")
    else
      tex.sprint("&")
    end
  end
  if tab[16] then
    tex.sprint("\\cellcolor[HTML]{C0C0C0}")
  end
end

return { ANDpocalipse = andpocalipse }