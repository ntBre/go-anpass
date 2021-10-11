anpass: *.go
	go build .

deploy: anpass
	scp -C anpass 'woods:Programs/anpass/.'
